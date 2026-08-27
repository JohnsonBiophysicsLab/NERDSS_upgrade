/*! \file test_check_for_structure_overlap.cpp
 *
 * ### Unit test for src/reactions/check_for_structure_overlap.cpp
 *
 * Function under test:
 *
 *     void check_for_structure_overlap(bool& cancelAssoc,
 *                                      const Complex& reactCom1,
 *                                      const Complex& reactCom2,
 *                                      const std::vector<Molecule>& moleculeList,
 *                                      const Parameters& params,
 *                                      const std::vector<MolTemplate>& molTemplateList)
 *
 * Behaviour, read directly from the implementation:
 *   - It computes `overlapTolerance = params.overlapSepLimit^2`.
 *   - It walks every member molecule of reactCom1 whose MolTemplate has
 *     `checkOverlap == true`, and compares it against every member molecule of
 *     reactCom2 whose MolTemplate also has `checkOverlap == true`.
 *   - The distance used is the *temporary association* coordinate
 *     (`Molecule::tmpComCoord`), NOT the permanent `comCoord`.
 *   - If the squared separation is **strictly less than** the tolerance, it
 *     sets `cancelAssoc = true` and returns immediately.
 *   - It never writes `false` into `cancelAssoc`; an incoming `true` survives.
 *
 * Every assertion below is written against exactly those semantics.
 */

#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (prefixed `cfso_` so they cannot collide with other test files)
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully-initialised Molecule for the overlap check.
 *
 * The function under test reads only `molTypeIndex` and `tmpComCoord`, but we
 * fill in `index` and `comCoord` too so the object is never left in a partially
 * initialised state.  `comCoord` is deliberately given a *different* value from
 * `tmpComCoord` in some tests to prove the routine really uses `tmpComCoord`.
 *
 * \param[in] index        index of this Molecule inside moleculeList
 * \param[in] molTypeIndex index into molTemplateList
 * \param[in] tmpCom       temporary association centre-of-mass coordinate
 * \param[in] permCom      permanent centre-of-mass coordinate
 */
Molecule cfso_make_molecule(int index, int molTypeIndex, const Coord& tmpCom,
                            const Coord& permCom)
{
    Molecule mol;
    mol.index = index;
    mol.molTypeIndex = molTypeIndex;
    mol.myComIndex = 0;
    mol.mass = 1.0;
    mol.isEmpty = false;
    mol.comCoord = permCom;
    mol.tmpComCoord = tmpCom;
    // No interfaces are needed: the routine only looks at the COM coordinates.
    mol.interfaceList.clear();
    mol.tmpICoords.clear();
    return mol;
}

/*! \brief Convenience overload: temporary and permanent coordinates identical. */
Molecule cfso_make_molecule(int index, int molTypeIndex, const Coord& com)
{
    return cfso_make_molecule(index, molTypeIndex, com, com);
}

/*! \brief Build a minimal MolTemplate carrying only the checkOverlap flag. */
MolTemplate cfso_make_template(int molTypeIndex, bool checkOverlap,
                               const std::string& name)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.checkOverlap = checkOverlap;
    temp.molName = name;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.D = Coord{ 1.0, 1.0, 1.0 };
    temp.Dr = Coord{ 0.01, 0.01, 0.01 };
    temp.interfaceList.clear();
    return temp;
}

/*! \brief Build a Complex that owns the supplied member-molecule indices. */
Complex cfso_make_complex(int index, const std::vector<int>& members)
{
    Complex com;
    com.index = index;
    com.memberList = members;
    com.comCoord = Coord{ 0.0, 0.0, 0.0 };
    com.D = Coord{ 1.0, 1.0, 1.0 };
    com.Dr = Coord{ 0.01, 0.01, 0.01 };
    com.mass = static_cast<double>(members.size());
    com.radius = 1.0;
    com.isEmpty = false;
    return com;
}

/*! \brief Parameters object holding just the overlap separation limit. */
Parameters cfso_make_params(double overlapSepLimit)
{
    Parameters params;
    params.overlapSepLimit = overlapSepLimit;
    params.timeStep = 1.0;
    return params;
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: two overlapping molecules (separation < overlapSepLimit) -> cancel.
// -----------------------------------------------------------------------------
void test_cfso_detects_overlap()
{
    std::cerr << "\n[TEST] test_cfso_detects_overlap\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    one molecule in each complex, both flagged\n"
              << "               checkOverlap=true, separated by 2.0 nm while\n"
              << "               overlapSepLimit is 10.0 nm.\n"
              << "  Pass:        cancelAssoc becomes true.\n";

    // Both molecule types must be checked for overlap.
    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B")
    };

    // Molecule 0 at the origin, molecule 1 only 2 nm away along x.
    std::vector<Molecule> moleculeList{
        cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
        cfso_make_molecule(1, 1, Coord{ 2.0, 0.0, 0.0 })
    };

    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });

    // Tolerance of 10 nm -> squared tolerance 100; the pair is at r^2 = 4.
    Parameters params = cfso_make_params(10.0);

    bool cancelAssoc{ false };
    std::cerr << "  Separation = 2.0 nm, tolerance = " << params.overlapSepLimit
              << " nm -> expect cancellation.\n";
    check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2, moleculeList,
                                params, molTemplateList);

    EXPECT_TRUE(cancelAssoc)
        << "Molecules 2 nm apart with a 10 nm overlap limit must cancel association";
    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc << '\n';
}

// -----------------------------------------------------------------------------
// Test 2: molecules further apart than the limit -> flag left untouched.
// -----------------------------------------------------------------------------
void test_cfso_no_overlap_when_far_apart()
{
    std::cerr << "\n[TEST] test_cfso_no_overlap_when_far_apart\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    the two molecules sit 50 nm apart while the\n"
              << "               overlap limit is only 10 nm.\n"
              << "  Pass:        cancelAssoc remains false.\n";

    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B")
    };

    std::vector<Molecule> moleculeList{
        cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
        cfso_make_molecule(1, 1, Coord{ 50.0, 0.0, 0.0 })
    };

    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });
    Parameters params = cfso_make_params(10.0);

    bool cancelAssoc{ false };
    std::cerr << "  Separation = 50.0 nm, tolerance = " << params.overlapSepLimit
              << " nm -> expect no cancellation.\n";
    check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2, moleculeList,
                                params, molTemplateList);

    EXPECT_FALSE(cancelAssoc)
        << "Well-separated molecules must not trigger an overlap cancellation";
    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc << '\n';
}

// -----------------------------------------------------------------------------
// Test 3: checkOverlap == false on either template disables the test entirely.
// -----------------------------------------------------------------------------
void test_cfso_respects_checkOverlap_flag()
{
    std::cerr << "\n[TEST] test_cfso_respects_checkOverlap_flag\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    molecules are coincident (r = 0) but at least\n"
              << "               one MolTemplate has checkOverlap = false.\n"
              << "  Pass:        cancelAssoc stays false in all three\n"
              << "               'flag disabled' permutations.\n";

    // Molecules are on top of each other -- the only thing that can stop a
    // cancellation here is the checkOverlap flag.
    std::vector<Molecule> moleculeList{
        cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
        cfso_make_molecule(1, 1, Coord{ 0.0, 0.0, 0.0 })
    };

    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });
    Parameters params = cfso_make_params(10.0);

    // Case (a): the first complex's molecule is not checked.
    {
        std::vector<MolTemplate> molTemplateList{
            cfso_make_template(0, false, "A"),
            cfso_make_template(1, true, "B")
        };
        bool cancelAssoc{ false };
        std::cerr << "  (a) reactCom1 template checkOverlap = false\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "Outer-loop molecule with checkOverlap=false must be skipped";
    }

    // Case (b): the second complex's molecule is not checked.
    {
        std::vector<MolTemplate> molTemplateList{
            cfso_make_template(0, true, "A"),
            cfso_make_template(1, false, "B")
        };
        bool cancelAssoc{ false };
        std::cerr << "  (b) reactCom2 template checkOverlap = false\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "Inner-loop molecule with checkOverlap=false must be skipped";
    }

    // Case (c): neither molecule is checked.
    {
        std::vector<MolTemplate> molTemplateList{
            cfso_make_template(0, false, "A"),
            cfso_make_template(1, false, "B")
        };
        bool cancelAssoc{ false };
        std::cerr << "  (c) both templates checkOverlap = false\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "With both flags off no pair may ever be compared";
    }
}

// -----------------------------------------------------------------------------
// Test 4: the comparison is strict (`r2 < tolerance`), so a separation exactly
//         equal to overlapSepLimit must NOT cancel.
// -----------------------------------------------------------------------------
void test_cfso_boundary_is_strictly_less_than()
{
    std::cerr << "\n[TEST] test_cfso_boundary_is_strictly_less_than\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    separation is exactly the overlap limit, then\n"
              << "               marginally smaller than the limit.\n"
              << "  Pass:        exactly-at-limit does not cancel (code uses\n"
              << "               `r2 < overlapTolerance`), just-inside does.\n";

    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B")
    };

    const double limit = 10.0;
    Parameters params = cfso_make_params(limit);

    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });

    // Case (a): r == limit exactly -> r2 == tolerance -> NOT strictly less.
    {
        std::vector<Molecule> moleculeList{
            cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
            cfso_make_molecule(1, 1, Coord{ limit, 0.0, 0.0 })
        };
        bool cancelAssoc{ false };
        std::cerr << "  (a) separation = " << limit << " nm (exactly the limit)\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "r2 == overlapTolerance is not '<' the tolerance, so no cancel";
    }

    // Case (b): a hair inside the limit -> must cancel.
    {
        const double justInside = limit - 1.0e-3;
        std::vector<Molecule> moleculeList{
            cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
            cfso_make_molecule(1, 1, Coord{ justInside, 0.0, 0.0 })
        };
        bool cancelAssoc{ false };
        std::cerr << "  (b) separation = " << justInside << " nm (just inside)\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_TRUE(cancelAssoc)
            << "A separation slightly below the limit must cancel association";
    }
}

// -----------------------------------------------------------------------------
// Test 5: the routine reads tmpComCoord, not comCoord.
// -----------------------------------------------------------------------------
void test_cfso_uses_tmp_coords_only()
{
    std::cerr << "\n[TEST] test_cfso_uses_tmp_coords_only\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    (a) permanent coords overlap but temporary\n"
              << "               coords are far apart, and (b) the reverse.\n"
              << "  Pass:        only the temporary coordinates influence the\n"
              << "               result.\n";

    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B")
    };
    Parameters params = cfso_make_params(10.0);
    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });

    // Case (a): permanent coords coincide, temporary coords 40 nm apart.
    {
        std::vector<Molecule> moleculeList{
            cfso_make_molecule(0, 0, /*tmp*/ Coord{ 0.0, 0.0, 0.0 },
                               /*perm*/ Coord{ 0.0, 0.0, 0.0 }),
            cfso_make_molecule(1, 1, /*tmp*/ Coord{ 40.0, 0.0, 0.0 },
                               /*perm*/ Coord{ 0.0, 0.0, 0.0 })
        };
        bool cancelAssoc{ false };
        std::cerr << "  (a) comCoord overlap, tmpComCoord far apart\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "Overlapping permanent coordinates must be ignored by this routine";
    }

    // Case (b): permanent coords far apart, temporary coords coincide.
    {
        std::vector<Molecule> moleculeList{
            cfso_make_molecule(0, 0, /*tmp*/ Coord{ 0.0, 0.0, 0.0 },
                               /*perm*/ Coord{ 0.0, 0.0, 0.0 }),
            cfso_make_molecule(1, 1, /*tmp*/ Coord{ 0.0, 0.0, 0.0 },
                               /*perm*/ Coord{ 40.0, 0.0, 0.0 })
        };
        bool cancelAssoc{ false };
        std::cerr << "  (b) comCoord far apart, tmpComCoord overlapping\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_TRUE(cancelAssoc)
            << "Overlapping temporary coordinates must trigger cancellation";
    }
}

// -----------------------------------------------------------------------------
// Test 6: many members per complex; a single overlapping pair is enough.
// -----------------------------------------------------------------------------
void test_cfso_multi_member_complexes()
{
    std::cerr << "\n[TEST] test_cfso_multi_member_complexes\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    three molecules per complex; first all pairs\n"
              << "               are well separated, then one deep pair is\n"
              << "               moved close together.\n"
              << "  Pass:        no cancel in the first case, cancel in the\n"
              << "               second (the double loop covers every pair).\n";

    // Four distinct molecule types, all checked for overlap.
    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B"),
        cfso_make_template(2, true, "C"),
        cfso_make_template(3, true, "D")
    };
    Parameters params = cfso_make_params(5.0);

    // Complex 1 owns molecules 0,1,2 ; complex 2 owns molecules 3,4,5.
    Complex reactCom1 = cfso_make_complex(0, { 0, 1, 2 });
    Complex reactCom2 = cfso_make_complex(1, { 3, 4, 5 });

    // Case (a): the two clusters are 100 nm apart in x -> no cancel.
    {
        std::vector<Molecule> moleculeList{
            cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
            cfso_make_molecule(1, 1, Coord{ 0.0, 10.0, 0.0 }),
            cfso_make_molecule(2, 2, Coord{ 0.0, 20.0, 0.0 }),
            cfso_make_molecule(3, 3, Coord{ 100.0, 0.0, 0.0 }),
            cfso_make_molecule(4, 0, Coord{ 100.0, 10.0, 0.0 }),
            cfso_make_molecule(5, 1, Coord{ 100.0, 20.0, 0.0 })
        };
        bool cancelAssoc{ false };
        std::cerr << "  (a) all 9 cross pairs separated by >= 100 nm\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "No cross pair is within 5 nm, so no cancellation is expected";
    }

    // Case (b): the *last* member of each complex is brought close together.
    // This exercises the innermost iteration of the nested loops.
    {
        std::vector<Molecule> moleculeList{
            cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
            cfso_make_molecule(1, 1, Coord{ 0.0, 10.0, 0.0 }),
            cfso_make_molecule(2, 2, Coord{ 0.0, 20.0, 0.0 }),
            cfso_make_molecule(3, 3, Coord{ 100.0, 0.0, 0.0 }),
            cfso_make_molecule(4, 0, Coord{ 100.0, 10.0, 0.0 }),
            cfso_make_molecule(5, 1, Coord{ 1.0, 20.0, 0.0 }) // 1 nm from mol 2
        };
        bool cancelAssoc{ false };
        std::cerr << "  (b) molecules 2 and 5 are only 1 nm apart\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_TRUE(cancelAssoc)
            << "A single overlapping cross pair anywhere must cancel association";
    }
}

// -----------------------------------------------------------------------------
// Test 7: an incoming `true` flag is never cleared by this routine.
// -----------------------------------------------------------------------------
void test_cfso_never_clears_incoming_flag()
{
    std::cerr << "\n[TEST] test_cfso_never_clears_incoming_flag\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    cancelAssoc arrives as true while the two\n"
              << "               complexes are far apart.\n"
              << "  Pass:        the flag is still true afterwards -- the code\n"
              << "               only ever writes `true`.\n";

    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B")
    };
    std::vector<Molecule> moleculeList{
        cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
        cfso_make_molecule(1, 1, Coord{ 500.0, 0.0, 0.0 })
    };
    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });
    Parameters params = cfso_make_params(10.0);

    bool cancelAssoc{ true }; // some earlier check already cancelled
    check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2, moleculeList,
                                params, molTemplateList);

    EXPECT_TRUE(cancelAssoc)
        << "check_for_structure_overlap must never reset a previously set flag";
    std::cerr << "  cancelAssoc = " << std::boolalpha << cancelAssoc << '\n';
}

// -----------------------------------------------------------------------------
// Test 8: empty member lists exercise the loops zero times and are harmless.
// -----------------------------------------------------------------------------
void test_cfso_empty_member_lists()
{
    std::cerr << "\n[TEST] test_cfso_empty_member_lists\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    one or both complexes have an empty memberList.\n"
              << "  Pass:        the routine returns without touching the flag.\n";

    std::vector<MolTemplate> molTemplateList{ cfso_make_template(0, true, "A") };
    std::vector<Molecule> moleculeList{
        cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 })
    };
    Parameters params = cfso_make_params(10.0);

    // (a) outer complex empty -> outer loop body never runs.
    {
        Complex emptyCom = cfso_make_complex(0, {});
        Complex fullCom = cfso_make_complex(1, { 0 });
        bool cancelAssoc{ false };
        std::cerr << "  (a) reactCom1.memberList is empty\n";
        check_for_structure_overlap(cancelAssoc, emptyCom, fullCom, moleculeList,
                                    params, molTemplateList);
        EXPECT_FALSE(cancelAssoc) << "Empty outer complex cannot cause an overlap";
    }

    // (b) inner complex empty -> inner loop body never runs.
    {
        Complex fullCom = cfso_make_complex(0, { 0 });
        Complex emptyCom = cfso_make_complex(1, {});
        bool cancelAssoc{ false };
        std::cerr << "  (b) reactCom2.memberList is empty\n";
        check_for_structure_overlap(cancelAssoc, fullCom, emptyCom, moleculeList,
                                    params, molTemplateList);
        EXPECT_FALSE(cancelAssoc) << "Empty inner complex cannot cause an overlap";
    }

    // (c) both empty.
    {
        Complex emptyCom1 = cfso_make_complex(0, {});
        Complex emptyCom2 = cfso_make_complex(1, {});
        bool cancelAssoc{ false };
        std::cerr << "  (c) both member lists are empty\n";
        check_for_structure_overlap(cancelAssoc, emptyCom1, emptyCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc) << "Two empty complexes cannot overlap";
    }
}

// -----------------------------------------------------------------------------
// Test 9: the tolerance genuinely scales with params.overlapSepLimit.
// -----------------------------------------------------------------------------
void test_cfso_tolerance_scales_with_parameter()
{
    std::cerr << "\n[TEST] test_cfso_tolerance_scales_with_parameter\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    the same 3 nm separation is evaluated against\n"
              << "               overlapSepLimit = 1 nm and 6 nm.\n"
              << "  Pass:        no cancel with the small limit, cancel with\n"
              << "               the large one.\n";

    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B")
    };

    // Separation is sqrt(3^2) = 3 nm along the x axis.
    std::vector<Molecule> moleculeList{
        cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
        cfso_make_molecule(1, 1, Coord{ 3.0, 0.0, 0.0 })
    };
    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });

    // Small limit: 3 nm > 1 nm -> no cancellation.
    {
        Parameters params = cfso_make_params(1.0);
        bool cancelAssoc{ false };
        std::cerr << "  (a) overlapSepLimit = 1.0 nm, separation = 3.0 nm\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "3 nm separation is outside a 1 nm overlap limit";
    }

    // Large limit: 3 nm < 6 nm -> cancellation.
    {
        Parameters params = cfso_make_params(6.0);
        bool cancelAssoc{ false };
        std::cerr << "  (b) overlapSepLimit = 6.0 nm, separation = 3.0 nm\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_TRUE(cancelAssoc)
            << "3 nm separation is inside a 6 nm overlap limit";
    }
}

// -----------------------------------------------------------------------------
// Test 10: the separation is a genuine 3-D distance, not a per-axis test.
// -----------------------------------------------------------------------------
void test_cfso_uses_full_3d_distance()
{
    std::cerr << "\n[TEST] test_cfso_uses_full_3d_distance\n"
              << "  Source file: src/reactions/check_for_structure_overlap.cpp\n"
              << "  Function:    check_for_structure_overlap\n"
              << "  Scenario:    a (3,4,0) offset gives |r| = 5 exactly; the\n"
              << "               limit is set just above and just below 5.\n"
              << "  Pass:        cancellation follows the 3-D distance, even\n"
              << "               though each individual axis offset is < limit.\n";

    std::vector<MolTemplate> molTemplateList{
        cfso_make_template(0, true, "A"),
        cfso_make_template(1, true, "B")
    };

    // dx=3, dy=4, dz=0 -> r = 5.0 exactly.
    std::vector<Molecule> moleculeList{
        cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
        cfso_make_molecule(1, 1, Coord{ 3.0, 4.0, 0.0 })
    };
    Complex reactCom1 = cfso_make_complex(0, { 0 });
    Complex reactCom2 = cfso_make_complex(1, { 1 });

    // (a) limit 4.5: each axis offset (3 and 4) is below 4.5, but the true
    //     distance of 5.0 is not -> must NOT cancel.
    {
        Parameters params = cfso_make_params(4.5);
        bool cancelAssoc{ false };
        std::cerr << "  (a) overlapSepLimit = 4.5 nm, true distance = 5.0 nm\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_FALSE(cancelAssoc)
            << "The Euclidean distance (5.0) exceeds the 4.5 nm limit";
    }

    // (b) limit 5.5: 5.0 < 5.5 -> cancel.
    {
        Parameters params = cfso_make_params(5.5);
        bool cancelAssoc{ false };
        std::cerr << "  (b) overlapSepLimit = 5.5 nm, true distance = 5.0 nm\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    moleculeList, params, molTemplateList);
        EXPECT_TRUE(cancelAssoc)
            << "The Euclidean distance (5.0) is inside the 5.5 nm limit";
    }

    // (c) a pure z-axis offset also has to be picked up.
    {
        std::vector<Molecule> zMoleculeList{
            cfso_make_molecule(0, 0, Coord{ 0.0, 0.0, 0.0 }),
            cfso_make_molecule(1, 1, Coord{ 0.0, 0.0, 2.0 })
        };
        Parameters params = cfso_make_params(5.0);
        bool cancelAssoc{ false };
        std::cerr << "  (c) offset entirely along z (2.0 nm), limit 5.0 nm\n";
        check_for_structure_overlap(cancelAssoc, reactCom1, reactCom2,
                                    zMoleculeList, params, molTemplateList);
        EXPECT_TRUE(cancelAssoc)
            << "The z component must contribute to the separation test";
    }
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers -- each named helper is run inside its own TEST so that
// a failure in one scenario does not stop the others from executing.
// -----------------------------------------------------------------------------
TEST(CheckForStructureOverlap, DetectsOverlap) { test_cfso_detects_overlap(); }
TEST(CheckForStructureOverlap, NoOverlapWhenFarApart) { test_cfso_no_overlap_when_far_apart(); }
TEST(CheckForStructureOverlap, RespectsCheckOverlapFlag) { test_cfso_respects_checkOverlap_flag(); }
TEST(CheckForStructureOverlap, BoundaryIsStrictlyLessThan) { test_cfso_boundary_is_strictly_less_than(); }
TEST(CheckForStructureOverlap, UsesTmpCoordsOnly) { test_cfso_uses_tmp_coords_only(); }
TEST(CheckForStructureOverlap, MultiMemberComplexes) { test_cfso_multi_member_complexes(); }
TEST(CheckForStructureOverlap, NeverClearsIncomingFlag) { test_cfso_never_clears_incoming_flag(); }
TEST(CheckForStructureOverlap, EmptyMemberLists) { test_cfso_empty_member_lists(); }
TEST(CheckForStructureOverlap, ToleranceScalesWithParameter) { test_cfso_tolerance_scales_with_parameter(); }
TEST(CheckForStructureOverlap, UsesFull3DDistance) { test_cfso_uses_full_3d_distance(); }