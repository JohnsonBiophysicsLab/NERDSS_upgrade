/*! \file test_check_for_structure_overlap_system.cpp
 *
 * ### Unit test for src/reactions/check_for_structure_overlap_system.cpp
 *
 * The single function under test is:
 *
 * \code
 * void check_for_structure_overlap_system(bool& flag,
 *                                         const Complex& reactCom1,
 *                                         const Complex& reactCom2,
 *                                         std::vector<Molecule>& moleculeList,
 *                                         const Parameters& params,
 *                                         const std::vector<MolTemplate>& molTemplateList,
 *                                         const std::vector<Complex>& complexList,
 *                                         const std::vector<ForwardRxn>& forwardRxns,
 *                                         const std::vector<BackRxn>& backRxns);
 * \endcode
 *
 * It walks every Complex in the system (skipping empty complexes and the two
 * complexes that are currently associating) and, for every pair of
 * "checkOverlap" molecules, compares
 *
 *   - the *actual* center of mass (comCoord) of the system molecule, against
 *   - the *temporary/association* center of mass (tmpComCoord) of a molecule in
 *     one of the two associating complexes.
 *
 * If that separation is smaller than params.overlapSepLimit the association is
 * cancelled by setting `flag = true` and returning immediately.  If the COMs are
 * further apart than overlapSepLimit but closer than the sum of the two molecule
 * radii, it falls through to an interface-level overlap check
 * (measure_overlap_protein_interfaces for the reactCom1 members,
 * measure_overlap_free_protein_interfaces for the reactCom2 members).
 *
 * Note that the routine only ever *sets* the flag to true -- it never clears it.
 *
 * The tests below build tiny, fully-initialized systems (molecules, complexes,
 * templates) so that every branch above is exercised, and print what is being
 * tested and why an assertion should pass.
 */

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

// -----------------------------------------------------------------------------
// Small helpers used to build a minimal but *complete* system.  Every object is
// fully initialized: the routines under test index into interfaceList,
// tmpICoords, memberList and molTemplateList without bounds checking, so any
// missing field would be a crash rather than a failed expectation.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a MolTemplate with exactly one interface that owns one state.
 *
 * \param[in] molTypeIndex index of this template inside molTemplateList
 * \param[in] name         molecule name (only used for readable output)
 * \param[in] checkOverlap whether this molecule type participates in the check
 * \param[in] radius       bounding radius used for the "are the interfaces
 *                         possibly close?" pre-screen
 */
MolTemplate cfsos_make_template(int molTypeIndex, const std::string& name, bool checkOverlap, double radius)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.checkOverlap = checkOverlap;
    temp.radius = radius;
    temp.mass = 1.0;
    temp.copies = 1;
    temp.D = Coord { 1.0, 1.0, 1.0 };
    temp.Dr = Coord { 0.01, 0.01, 0.01 };

    // One interface, with a single (default) state whose absolute index is 0 and
    // which participates in no reactions.  This keeps the interface-level checks
    // well defined while leaving the reaction lists empty.
    Interface iface;
    iface.index = 0;
    iface.name = name + "_iface";
    iface.iCoord = Coord { 0.0, 0.0, 1.0 };
    iface.stateList.emplace_back(0); // Interface::State(int index)
    temp.interfaceList.push_back(iface);

    return temp;
}

/*! \brief Build a Molecule with one interface, matching cfsos_make_template().
 *
 * \param[in] index        index of this molecule inside moleculeList
 * \param[in] comIndex     index of the parent complex
 * \param[in] molTypeIndex index of its MolTemplate
 * \param[in] com          the "real" center of mass coordinate
 * \param[in] ifaceCoord   the "real" interface coordinate
 */
Molecule cfsos_make_molecule(
    int index, int comIndex, int molTypeIndex, const Coord& com, const Coord& ifaceCoord)
{
    Molecule mol;
    mol.index = index;
    mol.id = index;
    mol.myComIndex = comIndex;
    mol.complexId = comIndex;
    mol.molTypeIndex = molTypeIndex;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.isLipid = false;
    mol.isImplicitLipid = false;
    mol.trajStatus = TrajStatus::none;
    mol.comCoord = com;

    Molecule::Iface iface;
    iface.coord = ifaceCoord;
    iface.index = 0; // absolute state index (matches template state index 0)
    iface.relIndex = 0; // relative interface index
    iface.stateIndex = 0; // first (only) state of the interface
    iface.stateIden = '\0';
    iface.molTypeIndex = molTypeIndex;
    iface.isBound = false;
    mol.interfaceList.push_back(iface);

    // Default the association ("tmp") coordinates to the real ones; tests that
    // care will override them with cfsos_set_tmp_coords().
    mol.tmpComCoord = com;
    mol.tmpICoords.push_back(ifaceCoord);

    return mol;
}

/*! \brief Overwrite the association (temporary) coordinates of a molecule.
 *
 * The molecules belonging to reactCom1/reactCom2 are compared using these
 * temporary coordinates, since they represent the post-association placement.
 */
void cfsos_set_tmp_coords(Molecule& mol, const Coord& tmpCom, const Coord& tmpIface)
{
    mol.tmpComCoord = tmpCom;
    mol.tmpICoords.clear();
    mol.tmpICoords.push_back(tmpIface);
}

/*! \brief Build a Complex owning the listed member molecules. */
Complex cfsos_make_complex(int index, const std::vector<int>& members, const Coord& com, bool isEmpty = false)
{
    Complex com1;
    com1.index = index;
    com1.id = index;
    com1.comCoord = com;
    com1.memberList = members;
    com1.mass = static_cast<double>(members.size());
    com1.radius = 1.0;
    com1.isEmpty = isEmpty;
    com1.OnSurface = false;
    com1.D = Coord { 1.0, 1.0, 1.0 };
    com1.Dr = Coord { 0.01, 0.01, 0.01 };
    com1.trajStatus = TrajStatus::none;
    return com1;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: Nothing is close to anything -> the flag must remain false.
// -----------------------------------------------------------------------------
void test_cfsos_no_overlap_when_far_apart()
{
    std::cerr << "\n[TEST] test_cfsos_no_overlap_when_far_apart\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    a third complex sits 100 nm away from both associating\n"
              << "               complexes; radii are tiny so the interface pre-screen\n"
              << "               is never entered.\n"
              << "  Pass:        flag stays false (association is NOT cancelled).\n";

    Parameters params;
    params.overlapSepLimit = 1.0; // COM-COM separations below 1 nm cancel

    // One molecule type, overlap-checked, with a very small bounding radius.
    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 0.1) };

    // molecule 0 -> complex 0 (reactCom1), molecule 1 -> complex 1 (reactCom2),
    // molecule 2 -> complex 2 (a bystander far away).
    std::vector<Molecule> moleculeList {
        cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
        cfsos_make_molecule(1, 1, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 1.0 }),
        cfsos_make_molecule(2, 2, 0, Coord { 100.0, 100.0, 100.0 }, Coord { 100.0, 100.0, 101.0 })
    };
    // The two associating molecules keep their tmp coords near the origin.
    cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 1.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 10.0, 0.0, 0.0 }),
        cfsos_make_complex(2, { 2 }, Coord { 100.0, 100.0, 100.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_FALSE(flag) << "Well separated complexes must not cancel association";
}

// -----------------------------------------------------------------------------
// Test 2: A system molecule sits on top of a reactCom1 member -> cancel.
// -----------------------------------------------------------------------------
void test_cfsos_com_overlap_with_react1_cancels()
{
    std::cerr << "\n[TEST] test_cfsos_com_overlap_with_react1_cancels\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    a bystander molecule is 0.5 nm from the *temporary* COM\n"
              << "               of the reactCom1 member while overlapSepLimit = 1.0 nm.\n"
              << "  Pass:        flag becomes true (association cancelled).\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 0.1) };

    std::vector<Molecule> moleculeList {
        cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
        cfsos_make_molecule(1, 1, 0, Coord { 50.0, 0.0, 0.0 }, Coord { 50.0, 0.0, 1.0 }),
        // bystander only 0.5 nm from moleculeList[0].tmpComCoord
        cfsos_make_molecule(2, 2, 0, Coord { 0.5, 0.0, 0.0 }, Coord { 0.5, 0.0, 1.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 50.0, 0.0, 0.0 }, Coord { 50.0, 0.0, 1.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 50.0, 0.0, 0.0 }),
        cfsos_make_complex(2, { 2 }, Coord { 0.5, 0.0, 0.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  COM separation = 0.5 nm, overlapSepLimit = " << params.overlapSepLimit
              << " nm; flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_TRUE(flag) << "A COM separation below overlapSepLimit must cancel association";
}

// -----------------------------------------------------------------------------
// Test 3: The same, but the overlap is with a member of reactCom2.
// -----------------------------------------------------------------------------
void test_cfsos_com_overlap_with_react2_cancels()
{
    std::cerr << "\n[TEST] test_cfsos_com_overlap_with_react2_cancels\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    same as before but the clash is with the reactCom2\n"
              << "               member (exercises the second inner loop).\n"
              << "  Pass:        flag becomes true.\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 0.1) };

    std::vector<Molecule> moleculeList {
        cfsos_make_molecule(0, 0, 0, Coord { -50.0, 0.0, 0.0 }, Coord { -50.0, 0.0, 1.0 }),
        cfsos_make_molecule(1, 1, 0, Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 1.0 }),
        // bystander 0.25 nm from moleculeList[1].tmpComCoord
        cfsos_make_molecule(2, 2, 0, Coord { 10.25, 0.0, 0.0 }, Coord { 10.25, 0.0, 1.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { -50.0, 0.0, 0.0 }, Coord { -50.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 10.0, 0.0, 0.0 }, Coord { 10.0, 0.0, 1.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { -50.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 10.0, 0.0, 0.0 }),
        cfsos_make_complex(2, { 2 }, Coord { 10.25, 0.0, 0.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  COM separation = 0.25 nm; flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_TRUE(flag) << "Overlap with a reactCom2 member must also cancel association";
}

// -----------------------------------------------------------------------------
// Test 4: Molecule types flagged checkOverlap == false are never compared.
// -----------------------------------------------------------------------------
void test_cfsos_checkOverlap_false_is_ignored()
{
    std::cerr << "\n[TEST] test_cfsos_checkOverlap_false_is_ignored\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    two sub-cases with perfectly overlapping COMs, but one\n"
              << "               of the two molecule types has checkOverlap = false.\n"
              << "  Pass:        flag stays false in both sub-cases.\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    // ---- Sub-case (a): the *system* molecule opts out of overlap checking ----
    {
        std::cerr << "  -> sub-case (a): the bystander (system) molecule has checkOverlap = false\n";
        std::vector<MolTemplate> molTemplateList {
            cfsos_make_template(0, "A_checked", true, 0.1), // associating molecules
            cfsos_make_template(1, "B_ignored", false, 0.1) // bystander
        };

        std::vector<Molecule> moleculeList {
            cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
            cfsos_make_molecule(1, 1, 0, Coord { 50.0, 0.0, 0.0 }, Coord { 50.0, 0.0, 1.0 }),
            cfsos_make_molecule(2, 2, 1, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 })
        };
        cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
        cfsos_set_tmp_coords(moleculeList[1], Coord { 50.0, 0.0, 0.0 }, Coord { 50.0, 0.0, 1.0 });

        std::vector<Complex> complexList {
            cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
            cfsos_make_complex(1, { 1 }, Coord { 50.0, 0.0, 0.0 }),
            cfsos_make_complex(2, { 2 }, Coord { 0.0, 0.0, 0.0 })
        };

        bool flag { false };
        check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
            molTemplateList, complexList, forwardRxns, backRxns);

        std::cerr << "     flag after call = " << std::boolalpha << flag << '\n';
        EXPECT_FALSE(flag) << "A bystander with checkOverlap == false must be skipped entirely";
    }

    // ---- Sub-case (b): the *associating* molecules opt out -------------------
    {
        std::cerr << "  -> sub-case (b): the associating molecules have checkOverlap = false\n";
        std::vector<MolTemplate> molTemplateList {
            cfsos_make_template(0, "A_ignored", false, 0.1), // associating molecules
            cfsos_make_template(1, "B_checked", true, 0.1) // bystander
        };

        std::vector<Molecule> moleculeList {
            cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
            cfsos_make_molecule(1, 1, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
            cfsos_make_molecule(2, 2, 1, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 })
        };
        cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
        cfsos_set_tmp_coords(moleculeList[1], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });

        std::vector<Complex> complexList {
            cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
            cfsos_make_complex(1, { 1 }, Coord { 0.0, 0.0, 0.0 }),
            cfsos_make_complex(2, { 2 }, Coord { 0.0, 0.0, 0.0 })
        };

        bool flag { false };
        check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
            molTemplateList, complexList, forwardRxns, backRxns);

        std::cerr << "     flag after call = " << std::boolalpha << flag << '\n';
        EXPECT_FALSE(flag) << "Associating molecules with checkOverlap == false must be skipped";
    }
}

// -----------------------------------------------------------------------------
// Test 5: Complexes marked isEmpty are skipped even if they overlap.
// -----------------------------------------------------------------------------
void test_cfsos_empty_complex_is_skipped()
{
    std::cerr << "\n[TEST] test_cfsos_empty_complex_is_skipped\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    the clashing bystander complex is marked isEmpty = true\n"
              << "               (a destroyed/void slot in complexList).\n"
              << "  Pass:        flag stays false because empty complexes are ignored.\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 0.1) };

    std::vector<Molecule> moleculeList {
        cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
        cfsos_make_molecule(1, 1, 0, Coord { 50.0, 0.0, 0.0 }, Coord { 50.0, 0.0, 1.0 }),
        cfsos_make_molecule(2, 2, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 50.0, 0.0, 0.0 }, Coord { 50.0, 0.0, 1.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 50.0, 0.0, 0.0 }),
        cfsos_make_complex(2, { 2 }, Coord { 0.0, 0.0, 0.0 }, /*isEmpty=*/true)
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_FALSE(flag) << "Empty (void) complexes must never trigger an overlap cancellation";
}

// -----------------------------------------------------------------------------
// Test 6: The two associating complexes are never compared to themselves.
// -----------------------------------------------------------------------------
void test_cfsos_reacting_complexes_are_skipped()
{
    std::cerr << "\n[TEST] test_cfsos_reacting_complexes_are_skipped\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    complexList holds ONLY the two associating complexes and\n"
              << "               their molecules sit at exactly the same coordinate.\n"
              << "  Pass:        flag stays false: c == reactCom1.index / reactCom2.index\n"
              << "               are explicitly skipped (that pair is handled elsewhere by\n"
              << "               check_for_structure_overlap()).\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 0.1) };

    std::vector<Molecule> moleculeList {
        cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
        cfsos_make_molecule(1, 1, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 0.0, 0.0, 0.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_FALSE(flag) << "The two reacting complexes must not be compared against themselves";
}

// -----------------------------------------------------------------------------
// Test 7: The routine never clears an already-set flag.
// -----------------------------------------------------------------------------
void test_cfsos_flag_is_never_cleared()
{
    std::cerr << "\n[TEST] test_cfsos_flag_is_never_cleared\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    caller passes flag = true (a previous check already\n"
              << "               cancelled the association) with no geometric overlap.\n"
              << "  Pass:        flag is still true on return (the routine only ever sets\n"
              << "               the flag, it never resets it).\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 0.1) };

    std::vector<Molecule> moleculeList {
        cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
        cfsos_make_molecule(1, 1, 0, Coord { 20.0, 0.0, 0.0 }, Coord { 20.0, 0.0, 1.0 }),
        cfsos_make_molecule(2, 2, 0, Coord { 200.0, 0.0, 0.0 }, Coord { 200.0, 0.0, 1.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 20.0, 0.0, 0.0 }, Coord { 20.0, 0.0, 1.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 20.0, 0.0, 0.0 }),
        cfsos_make_complex(2, { 2 }, Coord { 200.0, 0.0, 0.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { true }; // pretend a previous test already cancelled association
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_TRUE(flag) << "An incoming true flag must be preserved";
}

// -----------------------------------------------------------------------------
// Test 8: COMs closer than the summed radii but not overlapping -> the
//         interface-level check (measure_overlap_protein_interfaces) is invoked
//         for the reactCom1 members and must NOT produce a false positive when
//         the interfaces are ~10 nm apart.
// -----------------------------------------------------------------------------
void test_cfsos_react1_interface_prescreen_no_false_positive()
{
    std::cerr << "\n[TEST] test_cfsos_react1_interface_prescreen_no_false_positive\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    bounding radii are 5 nm each (sum^2 = 100) and the COM\n"
              << "               separation is 3 nm, so 1 < r^2 = 9 < 100: the routine\n"
              << "               falls through to measure_overlap_protein_interfaces for\n"
              << "               the reactCom1 member.  The interfaces themselves are\n"
              << "               ~10.4 nm apart.\n"
              << "  Pass:        flag stays false (no spurious cancellation).\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    // radius 5.0 -> (5+5)^2 = 100 nm^2 pre-screen window
    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 5.0) };

    std::vector<Molecule> moleculeList {
        // reactCom1 member: tmp COM at the origin, tmp interface 5 nm below it
        cfsos_make_molecule(0, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, -5.0 }),
        // reactCom2 member: parked far away so only the reactCom1 branch is used
        cfsos_make_molecule(1, 1, 0, Coord { 200.0, 200.0, 200.0 }, Coord { 200.0, 200.0, 205.0 }),
        // bystander: COM 3 nm away, interface 5 nm above it
        cfsos_make_molecule(2, 2, 0, Coord { 3.0, 0.0, 0.0 }, Coord { 3.0, 0.0, 5.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, -5.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 200.0, 200.0, 200.0 }, Coord { 200.0, 200.0, 205.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { 0.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 200.0, 200.0, 200.0 }),
        cfsos_make_complex(2, { 2 }, Coord { 3.0, 0.0, 0.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    const double ifaceSep = std::sqrt(3.0 * 3.0 + 10.0 * 10.0);
    std::cerr << "  COM separation = 3 nm, interface separation = " << ifaceSep
              << " nm; flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_FALSE(flag) << "Distant interfaces inside the radius pre-screen must not cancel association";
}

// -----------------------------------------------------------------------------
// Test 9: Same pre-screen path, but for a reactCom2 member, which routes through
//         measure_overlap_free_protein_interfaces (it additionally consults the
//         templates/reaction lists; here the interfaces take part in no
//         reactions, so nothing may be reported).
// -----------------------------------------------------------------------------
void test_cfsos_react2_interface_prescreen_no_false_positive()
{
    std::cerr << "\n[TEST] test_cfsos_react2_interface_prescreen_no_false_positive\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    mirror of the previous test, but with the near-by molecule\n"
              << "               belonging to reactCom2, so the free-interface variant\n"
              << "               (measure_overlap_free_protein_interfaces) is used.  The\n"
              << "               reaction lists are empty and the interfaces are 10+ nm apart.\n"
              << "  Pass:        flag stays false.\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 5.0) };

    std::vector<Molecule> moleculeList {
        // reactCom1 member: parked far away
        cfsos_make_molecule(0, 0, 0, Coord { -200.0, -200.0, -200.0 }, Coord { -200.0, -200.0, -195.0 }),
        // reactCom2 member: tmp COM 3 nm from the bystander, interface 5 nm below
        cfsos_make_molecule(1, 1, 0, Coord { 3.0, 0.0, 0.0 }, Coord { 3.0, 0.0, -5.0 }),
        // bystander at the origin with its interface 5 nm above
        cfsos_make_molecule(2, 2, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 5.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { -200.0, -200.0, -200.0 }, Coord { -200.0, -200.0, -195.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 3.0, 0.0, 0.0 }, Coord { 3.0, 0.0, -5.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0 }, Coord { -200.0, -200.0, -200.0 }),
        cfsos_make_complex(1, { 1 }, Coord { 3.0, 0.0, 0.0 }),
        cfsos_make_complex(2, { 2 }, Coord { 0.0, 0.0, 0.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_FALSE(flag) << "Free interfaces that cannot react and are far apart must not cancel association";
}

// -----------------------------------------------------------------------------
// Test 10: Multi-member complexes: the clash involves the *second* member of a
//          multi-molecule bystander complex and the *second* member of
//          reactCom1, verifying the nested loops cover every member.
// -----------------------------------------------------------------------------
void test_cfsos_multi_member_complexes()
{
    std::cerr << "\n[TEST] test_cfsos_multi_member_complexes\n"
              << "  Source file: check_for_structure_overlap_system.cpp\n"
              << "  Function:    check_for_structure_overlap_system\n"
              << "  Scenario:    reactCom1 holds two molecules and the bystander complex\n"
              << "               holds two molecules; only the second member of each pair\n"
              << "               clashes (0.1 nm apart).\n"
              << "  Pass:        flag becomes true, i.e. all member pairs are inspected.\n";

    Parameters params;
    params.overlapSepLimit = 1.0;

    std::vector<MolTemplate> molTemplateList { cfsos_make_template(0, "A", true, 0.1) };

    std::vector<Molecule> moleculeList {
        // complex 0 (reactCom1): molecules 0 and 1
        cfsos_make_molecule(0, 0, 0, Coord { -30.0, 0.0, 0.0 }, Coord { -30.0, 0.0, 1.0 }),
        cfsos_make_molecule(1, 0, 0, Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 }),
        // complex 1 (reactCom2): molecule 2, far from everything
        cfsos_make_molecule(2, 1, 0, Coord { 80.0, 0.0, 0.0 }, Coord { 80.0, 0.0, 1.0 }),
        // complex 2 (bystander): molecules 3 (far) and 4 (clashing with molecule 1)
        cfsos_make_molecule(3, 2, 0, Coord { 60.0, 60.0, 60.0 }, Coord { 60.0, 60.0, 61.0 }),
        cfsos_make_molecule(4, 2, 0, Coord { 0.1, 0.0, 0.0 }, Coord { 0.1, 0.0, 1.0 })
    };
    cfsos_set_tmp_coords(moleculeList[0], Coord { -30.0, 0.0, 0.0 }, Coord { -30.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[1], Coord { 0.0, 0.0, 0.0 }, Coord { 0.0, 0.0, 1.0 });
    cfsos_set_tmp_coords(moleculeList[2], Coord { 80.0, 0.0, 0.0 }, Coord { 80.0, 0.0, 1.0 });

    std::vector<Complex> complexList {
        cfsos_make_complex(0, { 0, 1 }, Coord { -15.0, 0.0, 0.0 }),
        cfsos_make_complex(1, { 2 }, Coord { 80.0, 0.0, 0.0 }),
        cfsos_make_complex(2, { 3, 4 }, Coord { 30.0, 30.0, 30.0 })
    };

    std::vector<ForwardRxn> forwardRxns {};
    std::vector<BackRxn> backRxns {};

    bool flag { false };
    check_for_structure_overlap_system(flag, complexList[0], complexList[1], moleculeList, params,
        molTemplateList, complexList, forwardRxns, backRxns);

    std::cerr << "  clashing pair separation = 0.1 nm; flag after call = " << std::boolalpha << flag << '\n';
    EXPECT_TRUE(flag) << "Every member of every complex must be inspected for overlap";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* function runs inside its own TEST so
// a failure in one does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(CheckForStructureOverlapSystem, NoOverlapWhenFarApart) { test_cfsos_no_overlap_when_far_apart(); }
TEST(CheckForStructureOverlapSystem, ComOverlapWithReact1Cancels) { test_cfsos_com_overlap_with_react1_cancels(); }
TEST(CheckForStructureOverlapSystem, ComOverlapWithReact2Cancels) { test_cfsos_com_overlap_with_react2_cancels(); }
TEST(CheckForStructureOverlapSystem, CheckOverlapFalseIsIgnored) { test_cfsos_checkOverlap_false_is_ignored(); }
TEST(CheckForStructureOverlapSystem, EmptyComplexIsSkipped) { test_cfsos_empty_complex_is_skipped(); }
TEST(CheckForStructureOverlapSystem, ReactingComplexesAreSkipped) { test_cfsos_reacting_complexes_are_skipped(); }
TEST(CheckForStructureOverlapSystem, FlagIsNeverCleared) { test_cfsos_flag_is_never_cleared(); }
TEST(CheckForStructureOverlapSystem, React1InterfacePrescreen) { test_cfsos_react1_interface_prescreen_no_false_positive(); }
TEST(CheckForStructureOverlapSystem, React2InterfacePrescreen) { test_cfsos_react2_interface_prescreen_no_false_positive(); }
TEST(CheckForStructureOverlapSystem, MultiMemberComplexes) { test_cfsos_multi_member_complexes(); }