/*! \file test_determine_shape_molecule.cpp
 *
 * ### Unit test for ../src/system_setup/determine_shape_molecule.cpp
 *
 * The function under test is:
 *
 *     void determine_shape_molecule(std::vector<MolTemplate>& molTemplateList)
 *
 * For every MolTemplate in the list it performs two distinct jobs:
 *
 *  1. Diffusion-constant sanitation
 *     - a lipid (isLipid) or implicit lipid (isImplicitLipid) is forced to
 *       have D.z == 0 (it lives on the 2D membrane),
 *     - a promoter (isPromoter) is forced to have D.y == 0 and D.z == 0
 *       (it lives on a 1D fiber). D.x is never touched.
 *
 *  2. Shape classification
 *     - isPoint is true when *every* interface coordinate equals the COM
 *       coordinate. Note that Coord::operator== compares the values after
 *       roundv(), i.e. after truncation to 4 decimal places, so differences
 *       smaller than 5e-5 compare as equal.
 *     - a point molecule is never a rod (isRod = false).
 *     - a non-point molecule with exactly one interface is always a rod.
 *     - a non-point molecule with several interfaces is a rod only when every
 *       pair of interfaces is co-linear with the COM (is_co_linear()).
 *     - an empty interface list leaves isPoint = true / isRod = false because
 *       the "all interfaces on the COM" loop never runs.
 *
 * Every test below prints what it is doing and what makes it pass.
 */

#include "classes/class_Coord.hpp"
#include "classes/class_MolTemplate.hpp"
#include "system_setup/system_setup.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed dsm_ so they cannot clash with other test files).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully initialized MolTemplate for the tests.
 *
 * Every member the function under test reads is set explicitly so that no
 * uninitialized data is ever touched.
 *
 * \param[in] name    name used in the warning messages.
 * \param[in] com     center-of-mass coordinate of the template.
 * \param[in] ifaces  list of interface coordinates to attach.
 * \return A MolTemplate ready to be handed to determine_shape_molecule().
 */
MolTemplate dsm_make_template(const std::string& name, const Coord& com,
    const std::vector<Coord>& ifaces)
{
    MolTemplate temp{};
    temp.molName = name;
    temp.comCoord = com;

    // Flags default to false; tests flip the ones they care about.
    temp.isLipid = false;
    temp.isImplicitLipid = false;
    temp.isPromoter = false;

    // Give the template well defined, non-zero diffusion constants so that a
    // change made by the function is unambiguous.
    temp.D = Coord{ 1.0, 2.0, 3.0 };
    temp.Dr = Coord{ 0.1, 0.2, 0.3 };

    // Deliberately set the shape flags to the "wrong" values so we can see
    // that the function really assigns them.
    temp.isPoint = false;
    temp.isRod = false;

    int index{ 0 };
    for (const auto& crd : ifaces) {
        temp.interfaceList.emplace_back("iface" + std::to_string(index), crd);
        temp.interfaceList.back().index = index;
        ++index;
    }

    return temp;
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: a lipid must have its z diffusion constant zeroed out.
// -----------------------------------------------------------------------------
void test_dsm_lipid_dz_is_zeroed()
{
    std::cerr << "\n[TEST] test_dsm_lipid_dz_is_zeroed\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    template flagged isLipid with D = (1,2,3).\n"
              << "  Pass:        D.z becomes exactly 0 while D.x/D.y are kept.\n";

    // A lipid whose single interface sits away from the COM.
    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("lipid", Coord{ 0.0, 0.0, 0.0 }, { Coord{ 0.0, 0.0, 1.0 } })
    };
    molTemplateList[0].isLipid = true;

    determine_shape_molecule(molTemplateList);

    EXPECT_DOUBLE_EQ(molTemplateList[0].D.z, 0.0)
        << "D.z of a lipid must be forced to zero";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.x, 1.0)
        << "D.x of a lipid must be left untouched";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.y, 2.0)
        << "D.y of a lipid must be left untouched";

    std::cerr << "  Resulting D = (" << molTemplateList[0].D.x << ", "
              << molTemplateList[0].D.y << ", " << molTemplateList[0].D.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 2: an implicit lipid must also have its z diffusion constant zeroed.
// -----------------------------------------------------------------------------
void test_dsm_implicit_lipid_dz_is_zeroed()
{
    std::cerr << "\n[TEST] test_dsm_implicit_lipid_dz_is_zeroed\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    template flagged isImplicitLipid (but NOT isLipid)\n"
              << "               with D = (1,2,3).\n"
              << "  Pass:        D.z becomes exactly 0, D.x/D.y unchanged.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("implicitLipid", Coord{ 0.0, 0.0, 0.0 }, { Coord{ 0.0, 0.0, 1.0 } })
    };
    molTemplateList[0].isImplicitLipid = true;

    determine_shape_molecule(molTemplateList);

    EXPECT_DOUBLE_EQ(molTemplateList[0].D.z, 0.0)
        << "D.z of an implicit lipid must be forced to zero";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.x, 1.0)
        << "D.x of an implicit lipid must be left untouched";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.y, 2.0)
        << "D.y of an implicit lipid must be left untouched";

    std::cerr << "  Resulting D = (" << molTemplateList[0].D.x << ", "
              << molTemplateList[0].D.y << ", " << molTemplateList[0].D.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 3: a promoter is confined to 1D, so D.y and D.z must both be zeroed.
// -----------------------------------------------------------------------------
void test_dsm_promoter_dy_dz_are_zeroed()
{
    std::cerr << "\n[TEST] test_dsm_promoter_dy_dz_are_zeroed\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    template flagged isPromoter with D = (1,2,3).\n"
              << "  Pass:        D.y and D.z become 0, D.x stays 1 (1D motion).\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("promoter", Coord{ 0.0, 0.0, 0.0 }, { Coord{ 1.0, 0.0, 0.0 } })
    };
    molTemplateList[0].isPromoter = true;

    determine_shape_molecule(molTemplateList);

    EXPECT_DOUBLE_EQ(molTemplateList[0].D.y, 0.0)
        << "D.y of a promoter must be forced to zero";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.z, 0.0)
        << "D.z of a promoter must be forced to zero";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.x, 1.0)
        << "D.x of a promoter must be left untouched (motion along the fiber)";

    std::cerr << "  Resulting D = (" << molTemplateList[0].D.x << ", "
              << molTemplateList[0].D.y << ", " << molTemplateList[0].D.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 4: an ordinary (non lipid / non promoter) molecule keeps all of its
//         diffusion constants.
// -----------------------------------------------------------------------------
void test_dsm_normal_molecule_keeps_diffusion()
{
    std::cerr << "\n[TEST] test_dsm_normal_molecule_keeps_diffusion\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    plain solution molecule, no special flags set.\n"
              << "  Pass:        D is returned untouched as (1,2,3).\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("solutionMol", Coord{ 0.0, 0.0, 0.0 }, { Coord{ 1.0, 1.0, 1.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_DOUBLE_EQ(molTemplateList[0].D.x, 1.0) << "D.x must not change";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.y, 2.0) << "D.y must not change";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.z, 3.0) << "D.z must not change";

    std::cerr << "  Resulting D = (" << molTemplateList[0].D.x << ", "
              << molTemplateList[0].D.y << ", " << molTemplateList[0].D.z << ")\n";
}

// -----------------------------------------------------------------------------
// Test 5: every interface sitting on the COM makes the molecule a point.
// -----------------------------------------------------------------------------
void test_dsm_point_molecule()
{
    std::cerr << "\n[TEST] test_dsm_point_molecule\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    two interfaces, both exactly on the COM.\n"
              << "  Pass:        isPoint == true and isRod == false.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("pointMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 0.0, 0.0, 0.0 }, Coord{ 0.0, 0.0, 0.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_TRUE(molTemplateList[0].isPoint)
        << "all interfaces coincide with the COM => isPoint";
    EXPECT_FALSE(molTemplateList[0].isRod)
        << "a point molecule is explicitly not a rod";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 6: a COM offset smaller than the roundv() resolution still counts as a
//         point, because Coord::operator== compares rounded (4 dp) values.
// -----------------------------------------------------------------------------
void test_dsm_point_within_rounding_tolerance()
{
    std::cerr << "\n[TEST] test_dsm_point_within_rounding_tolerance\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    interface offset by 1e-5 from the COM. Coord's\n"
              << "               operator== rounds to 4 decimal places (roundv),\n"
              << "               so 1e-5 rounds to 0 and the two compare equal.\n"
              << "  Pass:        isPoint == true, isRod == false.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("nearPointMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 1e-5, 0.0, 0.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_TRUE(molTemplateList[0].isPoint)
        << "sub-1e-4 offsets round to the COM and are treated as a point";
    EXPECT_FALSE(molTemplateList[0].isRod)
        << "a point molecule is explicitly not a rod";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 7: exactly one interface, off the COM, is always a rod.
// -----------------------------------------------------------------------------
void test_dsm_single_offset_interface_is_rod()
{
    std::cerr << "\n[TEST] test_dsm_single_offset_interface_is_rod\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    exactly one interface, displaced from the COM.\n"
              << "  Pass:        isPoint == false and isRod == true.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("singleIfaceMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 0.0, 0.0, 5.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_FALSE(molTemplateList[0].isPoint)
        << "an interface away from the COM means the molecule is not a point";
    EXPECT_TRUE(molTemplateList[0].isRod)
        << "a single off-COM interface is defined as a rod";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 8: two interfaces co-linear with the COM form a rod.
// -----------------------------------------------------------------------------
void test_dsm_two_colinear_interfaces_is_rod()
{
    std::cerr << "\n[TEST] test_dsm_two_colinear_interfaces_is_rod\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    COM at origin, interfaces at (1,0,0) and (-1,0,0);\n"
              << "               Heron's-formula area of that triple is exactly 0.\n"
              << "  Pass:        isPoint == false and isRod == true.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("rodMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 1.0, 0.0, 0.0 }, Coord{ -1.0, 0.0, 0.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_FALSE(molTemplateList[0].isPoint)
        << "interfaces are off the COM so it cannot be a point";
    EXPECT_TRUE(molTemplateList[0].isRod)
        << "COM and both interfaces are co-linear => rod";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 9: three interfaces that all lie on the same line are still a rod.
// -----------------------------------------------------------------------------
void test_dsm_three_colinear_interfaces_is_rod()
{
    std::cerr << "\n[TEST] test_dsm_three_colinear_interfaces_is_rod\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    COM at origin, interfaces at (0,0,1), (0,0,2)\n"
              << "               and (0,0,-3): every pair is co-linear with COM.\n"
              << "  Pass:        isPoint == false and isRod == true.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("longRodMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 0.0, 0.0, 1.0 }, Coord{ 0.0, 0.0, 2.0 }, Coord{ 0.0, 0.0, -3.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_FALSE(molTemplateList[0].isPoint) << "interfaces are off the COM";
    EXPECT_TRUE(molTemplateList[0].isRod)
        << "all interface pairs are co-linear with the COM => rod";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 10: two interfaces forming a triangle with the COM are not a rod.
// -----------------------------------------------------------------------------
void test_dsm_non_colinear_interfaces_is_not_rod()
{
    std::cerr << "\n[TEST] test_dsm_non_colinear_interfaces_is_not_rod\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    COM at origin with interfaces at (1,0,0) and\n"
              << "               (0,1,0): the triple encloses a finite area.\n"
              << "  Pass:        isPoint == false and isRod == false.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("bentMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 1.0, 0.0, 0.0 }, Coord{ 0.0, 1.0, 0.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_FALSE(molTemplateList[0].isPoint) << "interfaces are off the COM";
    EXPECT_FALSE(molTemplateList[0].isRod)
        << "a non co-linear interface pair rules out a rod";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 11: a bent molecule whose *first* pair is co-linear but whose later pair
//          is not must still end up as "not a rod".
// -----------------------------------------------------------------------------
void test_dsm_partially_colinear_is_not_rod()
{
    std::cerr << "\n[TEST] test_dsm_partially_colinear_is_not_rod\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    interfaces (0,0,1), (0,0,2) are co-linear with the\n"
              << "               COM, but the third interface (0,1,0) is not.\n"
              << "  Pass:        the pairwise loop must break out with isRod == false.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("partiallyStraightMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 0.0, 0.0, 1.0 }, Coord{ 0.0, 0.0, 2.0 }, Coord{ 0.0, 1.0, 0.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_FALSE(molTemplateList[0].isPoint) << "interfaces are off the COM";
    EXPECT_FALSE(molTemplateList[0].isRod)
        << "one non co-linear pair is enough to reject the rod classification";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 12: shape classification is done relative to the COM, not the origin.
// -----------------------------------------------------------------------------
void test_dsm_shape_relative_to_offset_com()
{
    std::cerr << "\n[TEST] test_dsm_shape_relative_to_offset_com\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    COM shifted to (5,5,5). Template A has both\n"
              << "               interfaces exactly on that COM; template B has\n"
              << "               interfaces symmetric about it along x.\n"
              << "  Pass:        A is a point (not rod), B is a rod (not point).\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("offsetPoint", Coord{ 5.0, 5.0, 5.0 },
            { Coord{ 5.0, 5.0, 5.0 }, Coord{ 5.0, 5.0, 5.0 } }),
        dsm_make_template("offsetRod", Coord{ 5.0, 5.0, 5.0 },
            { Coord{ 7.0, 5.0, 5.0 }, Coord{ 3.0, 5.0, 5.0 } })
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_TRUE(molTemplateList[0].isPoint)
        << "interfaces equal to the shifted COM => point";
    EXPECT_FALSE(molTemplateList[0].isRod) << "a point is never a rod";

    EXPECT_FALSE(molTemplateList[1].isPoint)
        << "interfaces displaced from the shifted COM => not a point";
    EXPECT_TRUE(molTemplateList[1].isRod)
        << "interfaces co-linear through the shifted COM => rod";

    std::cerr << "  offsetPoint: isPoint = " << std::boolalpha
              << molTemplateList[0].isPoint << ", isRod = " << molTemplateList[0].isRod
              << "\n  offsetRod:   isPoint = " << molTemplateList[1].isPoint
              << ", isRod = " << molTemplateList[1].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 13: a template with no interfaces at all.
// -----------------------------------------------------------------------------
void test_dsm_empty_interface_list()
{
    std::cerr << "\n[TEST] test_dsm_empty_interface_list\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    template carrying an empty interfaceList. The\n"
              << "               'all interfaces on the COM' loop never executes,\n"
              << "               so isPoint stays true and the function returns\n"
              << "               early with isRod = false.\n"
              << "  Pass:        isPoint == true and isRod == false (no crash).\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("noIfaceMol", Coord{ 0.0, 0.0, 0.0 }, {})
    };

    determine_shape_molecule(molTemplateList);

    EXPECT_TRUE(molTemplateList[0].isPoint)
        << "an empty interface list leaves isPoint at its initialized value";
    EXPECT_FALSE(molTemplateList[0].isRod)
        << "the point branch resets isRod to false";

    std::cerr << "  isPoint = " << std::boolalpha << molTemplateList[0].isPoint
              << ", isRod = " << molTemplateList[0].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 14: an empty template list is handled gracefully (nothing to do).
// -----------------------------------------------------------------------------
void test_dsm_empty_template_list()
{
    std::cerr << "\n[TEST] test_dsm_empty_template_list\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    an empty molTemplateList is passed in.\n"
              << "  Pass:        the call returns without touching anything.\n";

    std::vector<MolTemplate> molTemplateList{};

    determine_shape_molecule(molTemplateList);

    EXPECT_TRUE(molTemplateList.empty())
        << "the function must not add or remove templates";

    std::cerr << "  molTemplateList size after call = " << molTemplateList.size() << '\n';
}

// -----------------------------------------------------------------------------
// Test 15: every template in a mixed list is processed independently in one
//          single call (diffusion sanitation + shape classification together).
// -----------------------------------------------------------------------------
void test_dsm_processes_whole_list()
{
    std::cerr << "\n[TEST] test_dsm_processes_whole_list\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    one call with a lipid rod, a promoter point and\n"
              << "               a plain bent solution molecule.\n"
              << "  Pass:        each entry receives its own correct diffusion\n"
              << "               constants and its own shape flags.\n";

    std::vector<MolTemplate> molTemplateList{
        // [0] lipid, rod-shaped (two interfaces co-linear through the COM).
        dsm_make_template("lipidRod", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 2.0, 0.0, 0.0 }, Coord{ -4.0, 0.0, 0.0 } }),
        // [1] promoter whose single interface sits on the COM => point.
        dsm_make_template("promoterPoint", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 0.0, 0.0, 0.0 } }),
        // [2] plain molecule with a bent geometry => neither point nor rod.
        dsm_make_template("bentSolutionMol", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 1.0, 0.0, 0.0 }, Coord{ 0.0, 0.0, 1.0 } })
    };
    molTemplateList[0].isLipid = true;
    molTemplateList[1].isPromoter = true;

    determine_shape_molecule(molTemplateList);

    // --- entry [0]: lipid rod -------------------------------------------------
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.z, 0.0) << "lipid D.z zeroed";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.y, 2.0) << "lipid D.y untouched";
    EXPECT_FALSE(molTemplateList[0].isPoint) << "lipid interfaces are off the COM";
    EXPECT_TRUE(molTemplateList[0].isRod) << "lipid interfaces are co-linear";

    // --- entry [1]: promoter point -------------------------------------------
    EXPECT_DOUBLE_EQ(molTemplateList[1].D.y, 0.0) << "promoter D.y zeroed";
    EXPECT_DOUBLE_EQ(molTemplateList[1].D.z, 0.0) << "promoter D.z zeroed";
    EXPECT_DOUBLE_EQ(molTemplateList[1].D.x, 1.0) << "promoter D.x untouched";
    EXPECT_TRUE(molTemplateList[1].isPoint) << "promoter interface sits on the COM";
    EXPECT_FALSE(molTemplateList[1].isRod) << "a point is never a rod";

    // --- entry [2]: bent solution molecule ------------------------------------
    EXPECT_DOUBLE_EQ(molTemplateList[2].D.x, 1.0) << "plain molecule D.x untouched";
    EXPECT_DOUBLE_EQ(molTemplateList[2].D.y, 2.0) << "plain molecule D.y untouched";
    EXPECT_DOUBLE_EQ(molTemplateList[2].D.z, 3.0) << "plain molecule D.z untouched";
    EXPECT_FALSE(molTemplateList[2].isPoint) << "interfaces are off the COM";
    EXPECT_FALSE(molTemplateList[2].isRod) << "interfaces are not co-linear";

    std::cerr << "  [0] lipidRod        : isPoint = " << std::boolalpha
              << molTemplateList[0].isPoint << ", isRod = " << molTemplateList[0].isRod
              << ", D.z = " << molTemplateList[0].D.z << '\n'
              << "  [1] promoterPoint   : isPoint = " << molTemplateList[1].isPoint
              << ", isRod = " << molTemplateList[1].isRod
              << ", D = (" << molTemplateList[1].D.x << "," << molTemplateList[1].D.y
              << "," << molTemplateList[1].D.z << ")\n"
              << "  [2] bentSolutionMol : isPoint = " << molTemplateList[2].isPoint
              << ", isRod = " << molTemplateList[2].isRod << '\n';
}

// -----------------------------------------------------------------------------
// Test 16: calling the function twice on the same data must be idempotent.
// -----------------------------------------------------------------------------
void test_dsm_is_idempotent()
{
    std::cerr << "\n[TEST] test_dsm_is_idempotent\n"
              << "  Source file: system_setup/determine_shape_molecule.cpp\n"
              << "  Function:    determine_shape_molecule()\n"
              << "  Scenario:    the same lipid rod is classified twice in a row.\n"
              << "  Pass:        the second call reproduces the first result exactly.\n";

    std::vector<MolTemplate> molTemplateList{
        dsm_make_template("idempotentLipidRod", Coord{ 0.0, 0.0, 0.0 },
            { Coord{ 3.0, 0.0, 0.0 }, Coord{ -3.0, 0.0, 0.0 } })
    };
    molTemplateList[0].isLipid = true;

    determine_shape_molecule(molTemplateList);
    const bool firstIsPoint{ molTemplateList[0].isPoint };
    const bool firstIsRod{ molTemplateList[0].isRod };
    const double firstDz{ molTemplateList[0].D.z };

    determine_shape_molecule(molTemplateList);

    EXPECT_EQ(molTemplateList[0].isPoint, firstIsPoint)
        << "isPoint must be stable across repeated calls";
    EXPECT_EQ(molTemplateList[0].isRod, firstIsRod)
        << "isRod must be stable across repeated calls";
    EXPECT_DOUBLE_EQ(molTemplateList[0].D.z, firstDz)
        << "D.z must remain zero after the second call";

    std::cerr << "  after 2 calls: isPoint = " << std::boolalpha
              << molTemplateList[0].isPoint << ", isRod = " << molTemplateList[0].isRod
              << ", D.z = " << molTemplateList[0].D.z << '\n';
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each helper is invoked from its own TEST so that a
// failure in one scenario does not prevent the remaining ones from running.
// -----------------------------------------------------------------------------
TEST(DetermineShapeMolecule, LipidDzIsZeroed) { test_dsm_lipid_dz_is_zeroed(); }
TEST(DetermineShapeMolecule, ImplicitLipidDzIsZeroed) { test_dsm_implicit_lipid_dz_is_zeroed(); }
TEST(DetermineShapeMolecule, PromoterDyDzAreZeroed) { test_dsm_promoter_dy_dz_are_zeroed(); }
TEST(DetermineShapeMolecule, NormalMoleculeKeepsDiffusion) { test_dsm_normal_molecule_keeps_diffusion(); }
TEST(DetermineShapeMolecule, PointMolecule) { test_dsm_point_molecule(); }
TEST(DetermineShapeMolecule, PointWithinRoundingTolerance) { test_dsm_point_within_rounding_tolerance(); }
TEST(DetermineShapeMolecule, SingleOffsetInterfaceIsRod) { test_dsm_single_offset_interface_is_rod(); }
TEST(DetermineShapeMolecule, TwoColinearInterfacesIsRod) { test_dsm_two_colinear_interfaces_is_rod(); }
TEST(DetermineShapeMolecule, ThreeColinearInterfacesIsRod) { test_dsm_three_colinear_interfaces_is_rod(); }
TEST(DetermineShapeMolecule, NonColinearInterfacesIsNotRod) { test_dsm_non_colinear_interfaces_is_not_rod(); }
TEST(DetermineShapeMolecule, PartiallyColinearIsNotRod) { test_dsm_partially_colinear_is_not_rod(); }
TEST(DetermineShapeMolecule, ShapeRelativeToOffsetCom) { test_dsm_shape_relative_to_offset_com(); }
TEST(DetermineShapeMolecule, EmptyInterfaceList) { test_dsm_empty_interface_list(); }
TEST(DetermineShapeMolecule, EmptyTemplateList) { test_dsm_empty_template_list(); }
TEST(DetermineShapeMolecule, ProcessesWholeList) { test_dsm_processes_whole_list(); }
TEST(DetermineShapeMolecule, IsIdempotent) { test_dsm_is_idempotent(); }