/*! \file test_set_exclude_volume_bound.cpp
 *
 * ### Unit test for src/system_setup/set_exclude_volume_bound.cpp
 *
 * Function under test:
 *
 *     void set_exclude_volume_bound(std::vector<ForwardRxn>& forwardRxns,
 *                                   std::vector<MolTemplate>& molTemplateList);
 *
 * ## What the function actually does (read from the implementation)
 *
 * It loops over every ForwardRxn.  For each reaction whose `excludeVolumeBound`
 * flag is `true` it takes the *first two* entries of `reactantListNew` and
 * cross-registers them on the two participating MolTemplate interfaces:
 *
 *   - reactant-1's interface gets reactant-2's molTypeIndex appended to
 *     `excludeVolumeBoundList`, and vice-versa,
 *   - reactant-1's interface gets reactant-2's relIfaceIndex appended to
 *     `excludeVolumeBoundIfaceList`, and vice-versa,
 *   - *both* interfaces get the reaction's `bindRadius` appended to
 *     `excludeRadiusList`,
 *   - *both* interfaces get the reaction's `relRxnIndex` appended to
 *     `excludeVolumeBoundReactList`,
 *   - and both parent MolTemplates get `excludeVolumeBound = true`.
 *
 * Important behavioural details that the assertions below rely on:
 *   - The function **appends**; it never clears pre-existing entries.
 *   - Reactions with `excludeVolumeBound == false` are skipped entirely, so
 *     nothing at all is touched for them.
 *   - A self-reaction (same molTypeIndex *and* same relIfaceIndex on both
 *     sides) writes to the *same* interface twice, so each list gains two
 *     entries.
 *   - `absIfaceIndex` is read into local variables but is never used, so it
 *     cannot influence the result.
 *
 * All objects handed to the function are fully initialized here, because the
 * implementation indexes `molTemplateList[...]` and `.interfaceList[...]`
 * without any bounds checking.
 */

#include "classes/class_MolTemplate.hpp"
#include "classes/class_Rxns.hpp"
#include "system_setup/system_setup.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (uniquely prefixed with `sevb_` so they cannot collide with
// helpers defined by other test translation units in the same binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully initialized MolTemplate with `numIfaces` interfaces.
 *
 * Each interface gets a relative index, a name and one default state so that
 * the object looks like something the parser would have produced.  The four
 * exclude-volume vectors start out empty, which is what the function under
 * test expects to append to.
 *
 * \param[in] molTypeIndex index this template will occupy in molTemplateList
 * \param[in] name         molecule name (used to build interface names)
 * \param[in] numIfaces    how many interfaces to create
 */
MolTemplate sevb_make_mol_template(int molTypeIndex, const std::string& name,
                                   int numIfaces)
{
    MolTemplate temp;
    temp.molTypeIndex = molTypeIndex;
    temp.molName = name;
    temp.copies = 1;
    temp.mass = 1.0;
    temp.radius = 1.0;
    temp.excludeVolumeBound = false; // the function may flip this to true

    for (int i = 0; i < numIfaces; ++i) {
        Interface iface;
        iface.index = i;
        iface.name = name + "_iface" + std::to_string(i);
        // Give the interface a single default state, as the parser would.
        iface.stateList.emplace_back(iface.name, i);
        temp.interfaceList.push_back(iface);
    }
    return temp;
}

/*! \brief Build a minimal bimolecular ForwardRxn between two interfaces.
 *
 * Only the fields that `set_exclude_volume_bound` reads are meaningful:
 * `excludeVolumeBound`, `bindRadius`, `relRxnIndex` and the first two entries
 * of `reactantListNew` (their molTypeIndex / relIfaceIndex / absIfaceIndex).
 *
 * \param[in] excludeVolumeBound whether this reaction opts into the feature
 * \param[in] bindRadius         binding radius appended to excludeRadiusList
 * \param[in] relRxnIndex        index appended to excludeVolumeBoundReactList
 * \param[in] mol1,rel1,abs1     reactant 1 molecule / relative / absolute index
 * \param[in] mol2,rel2,abs2     reactant 2 molecule / relative / absolute index
 */
ForwardRxn sevb_make_rxn(bool excludeVolumeBound, double bindRadius,
                         int relRxnIndex, int mol1, int rel1, int abs1,
                         int mol2, int rel2, int abs2)
{
    ForwardRxn rxn;
    rxn.rxnType = ReactionType::bimolecular;
    rxn.excludeVolumeBound = excludeVolumeBound;
    rxn.bindRadius = bindRadius;
    rxn.relRxnIndex = relRxnIndex;
    rxn.absRxnIndex = relRxnIndex;

    // RxnIface(ifaceName, molTypeIndex, absIfaceIndex, relIfaceIndex,
    //          requiresState, requiresInteraction)
    rxn.reactantListNew.emplace_back("react1", mol1, abs1, rel1, '\0', false);
    rxn.reactantListNew.emplace_back("react2", mol2, abs2, rel2, '\0', false);

    return rxn;
}

/*! \brief Convenience printer used by the verbose console output. */
void sevb_dump_iface(const MolTemplate& temp, int ifaceIndex)
{
    const Interface& iface = temp.interfaceList[ifaceIndex];
    std::cerr << "      " << temp.molName << ".interfaceList[" << ifaceIndex
              << "]  molList={";
    for (size_t i = 0; i < iface.excludeVolumeBoundList.size(); ++i)
        std::cerr << (i ? "," : "") << iface.excludeVolumeBoundList[i];
    std::cerr << "}  ifaceList={";
    for (size_t i = 0; i < iface.excludeVolumeBoundIfaceList.size(); ++i)
        std::cerr << (i ? "," : "") << iface.excludeVolumeBoundIfaceList[i];
    std::cerr << "}  radii={";
    for (size_t i = 0; i < iface.excludeRadiusList.size(); ++i)
        std::cerr << (i ? "," : "") << iface.excludeRadiusList[i];
    std::cerr << "}  rxnList={";
    for (size_t i = 0; i < iface.excludeVolumeBoundReactList.size(); ++i)
        std::cerr << (i ? "," : "") << iface.excludeVolumeBoundReactList[i];
    std::cerr << "}\n";
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: an empty reaction list must leave every MolTemplate untouched.
// -----------------------------------------------------------------------------
void sevb_test_empty_reaction_list()
{
    std::cerr << "\n[TEST] sevb_test_empty_reaction_list\n"
              << "  Source file: src/system_setup/set_exclude_volume_bound.cpp\n"
              << "  Function:    set_exclude_volume_bound\n"
              << "  Scenario:    forwardRxns is empty.\n"
              << "  Criteria:    all four exclude-volume vectors stay empty and\n"
              << "               MolTemplate::excludeVolumeBound stays false.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(sevb_make_mol_template(0, "A", 2));
    molTemplateList.push_back(sevb_make_mol_template(1, "B", 2));

    std::vector<ForwardRxn> forwardRxns; // deliberately empty

    std::cerr << "  Calling set_exclude_volume_bound with 0 reactions...\n";
    set_exclude_volume_bound(forwardRxns, molTemplateList);

    for (const auto& temp : molTemplateList) {
        EXPECT_FALSE(temp.excludeVolumeBound)
            << "MolTemplate " << temp.molName
            << " must not be flagged when there are no reactions";
        for (const auto& iface : temp.interfaceList) {
            EXPECT_TRUE(iface.excludeVolumeBoundList.empty())
                << "excludeVolumeBoundList should still be empty";
            EXPECT_TRUE(iface.excludeVolumeBoundIfaceList.empty())
                << "excludeVolumeBoundIfaceList should still be empty";
            EXPECT_TRUE(iface.excludeRadiusList.empty())
                << "excludeRadiusList should still be empty";
            EXPECT_TRUE(iface.excludeVolumeBoundReactList.empty())
                << "excludeVolumeBoundReactList should still be empty";
        }
    }
    std::cerr << "  All templates verified untouched.\n";
}

// -----------------------------------------------------------------------------
// Test 2: reactions that do NOT opt in are skipped entirely.
// -----------------------------------------------------------------------------
void sevb_test_flag_false_is_skipped()
{
    std::cerr << "\n[TEST] sevb_test_flag_false_is_skipped\n"
              << "  Source file: src/system_setup/set_exclude_volume_bound.cpp\n"
              << "  Function:    set_exclude_volume_bound\n"
              << "  Scenario:    two reactions, both with excludeVolumeBound=false.\n"
              << "  Criteria:    nothing is appended and no template is flagged.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(sevb_make_mol_template(0, "A", 2));
    molTemplateList.push_back(sevb_make_mol_template(1, "B", 2));

    std::vector<ForwardRxn> forwardRxns;
    // A(iface0) + B(iface1), opted OUT
    forwardRxns.push_back(sevb_make_rxn(false, 1.0, 0, 0, 0, 10, 1, 1, 11));
    // A(iface1) + B(iface0), opted OUT
    forwardRxns.push_back(sevb_make_rxn(false, 2.0, 1, 0, 1, 12, 1, 0, 13));

    std::cerr << "  Calling set_exclude_volume_bound with 2 opted-out reactions...\n";
    set_exclude_volume_bound(forwardRxns, molTemplateList);

    for (const auto& temp : molTemplateList) {
        EXPECT_FALSE(temp.excludeVolumeBound)
            << "MolTemplate " << temp.molName
            << " must not be flagged by an opted-out reaction";
        for (const auto& iface : temp.interfaceList) {
            EXPECT_EQ(iface.excludeVolumeBoundList.size(), 0u)
                << "no partner molecule index should have been recorded";
            EXPECT_EQ(iface.excludeVolumeBoundIfaceList.size(), 0u)
                << "no partner interface index should have been recorded";
            EXPECT_EQ(iface.excludeRadiusList.size(), 0u)
                << "no bind radius should have been recorded";
            EXPECT_EQ(iface.excludeVolumeBoundReactList.size(), 0u)
                << "no reaction index should have been recorded";
        }
    }
    std::cerr << "  Confirmed the opted-out reactions produced no side effects.\n";
}

// -----------------------------------------------------------------------------
// Test 3: a single opted-in reaction cross-registers both partners.
// -----------------------------------------------------------------------------
void sevb_test_single_reaction_cross_registers()
{
    std::cerr << "\n[TEST] sevb_test_single_reaction_cross_registers\n"
              << "  Source file: src/system_setup/set_exclude_volume_bound.cpp\n"
              << "  Function:    set_exclude_volume_bound\n"
              << "  Scenario:    one reaction A(iface0) + B(iface1),\n"
              << "               excludeVolumeBound=true, bindRadius=3.5,\n"
              << "               relRxnIndex=7.\n"
              << "  Criteria:    A.iface0 records molecule 1 / iface 1 / 3.5 / 7,\n"
              << "               B.iface1 records molecule 0 / iface 0 / 3.5 / 7,\n"
              << "               both templates get excludeVolumeBound=true, and\n"
              << "               the uninvolved interfaces stay empty.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(sevb_make_mol_template(0, "A", 2));
    molTemplateList.push_back(sevb_make_mol_template(1, "B", 2));

    std::vector<ForwardRxn> forwardRxns;
    // mol 0, relIface 0, absIface 10  <->  mol 1, relIface 1, absIface 11
    forwardRxns.push_back(sevb_make_rxn(true, 3.5, 7, 0, 0, 10, 1, 1, 11));

    std::cerr << "  Calling set_exclude_volume_bound...\n";
    set_exclude_volume_bound(forwardRxns, molTemplateList);

    const Interface& aIface0 = molTemplateList[0].interfaceList[0];
    const Interface& bIface1 = molTemplateList[1].interfaceList[1];

    std::cerr << "  Resulting interface contents:\n";
    sevb_dump_iface(molTemplateList[0], 0);
    sevb_dump_iface(molTemplateList[1], 1);

    // --- Reactant 1 side: A.iface0 should point at molecule B, interface 1 ---
    ASSERT_EQ(aIface0.excludeVolumeBoundList.size(), 1u)
        << "A.iface0 should have exactly one partner molecule recorded";
    EXPECT_EQ(aIface0.excludeVolumeBoundList[0], 1)
        << "A.iface0 partner molecule index should be 1 (B)";

    ASSERT_EQ(aIface0.excludeVolumeBoundIfaceList.size(), 1u)
        << "A.iface0 should have exactly one partner interface recorded";
    EXPECT_EQ(aIface0.excludeVolumeBoundIfaceList[0], 1)
        << "A.iface0 partner interface index should be 1";

    ASSERT_EQ(aIface0.excludeRadiusList.size(), 1u)
        << "A.iface0 should have exactly one exclude radius recorded";
    EXPECT_DOUBLE_EQ(aIface0.excludeRadiusList[0], 3.5)
        << "A.iface0 exclude radius should be the reaction bindRadius (3.5)";

    ASSERT_EQ(aIface0.excludeVolumeBoundReactList.size(), 1u)
        << "A.iface0 should have exactly one reaction index recorded";
    EXPECT_EQ(aIface0.excludeVolumeBoundReactList[0], 7)
        << "A.iface0 reaction index should be relRxnIndex (7)";

    // --- Reactant 2 side: B.iface1 should point at molecule A, interface 0 ---
    ASSERT_EQ(bIface1.excludeVolumeBoundList.size(), 1u)
        << "B.iface1 should have exactly one partner molecule recorded";
    EXPECT_EQ(bIface1.excludeVolumeBoundList[0], 0)
        << "B.iface1 partner molecule index should be 0 (A)";

    ASSERT_EQ(bIface1.excludeVolumeBoundIfaceList.size(), 1u)
        << "B.iface1 should have exactly one partner interface recorded";
    EXPECT_EQ(bIface1.excludeVolumeBoundIfaceList[0], 0)
        << "B.iface1 partner interface index should be 0";

    ASSERT_EQ(bIface1.excludeRadiusList.size(), 1u)
        << "B.iface1 should have exactly one exclude radius recorded";
    EXPECT_DOUBLE_EQ(bIface1.excludeRadiusList[0], 3.5)
        << "B.iface1 exclude radius should be the reaction bindRadius (3.5)";

    ASSERT_EQ(bIface1.excludeVolumeBoundReactList.size(), 1u)
        << "B.iface1 should have exactly one reaction index recorded";
    EXPECT_EQ(bIface1.excludeVolumeBoundReactList[0], 7)
        << "B.iface1 reaction index should be relRxnIndex (7)";

    // --- Template-level flags ---
    EXPECT_TRUE(molTemplateList[0].excludeVolumeBound)
        << "MolTemplate A should be flagged excludeVolumeBound";
    EXPECT_TRUE(molTemplateList[1].excludeVolumeBound)
        << "MolTemplate B should be flagged excludeVolumeBound";

    // --- Interfaces not involved in the reaction must remain untouched ---
    EXPECT_TRUE(molTemplateList[0].interfaceList[1].excludeVolumeBoundList.empty())
        << "A.iface1 was not part of the reaction and must stay empty";
    EXPECT_TRUE(molTemplateList[1].interfaceList[0].excludeVolumeBoundList.empty())
        << "B.iface0 was not part of the reaction and must stay empty";
}

// -----------------------------------------------------------------------------
// Test 4: several reactions accumulate on the same interface, in order, and
//         opted-out reactions interleaved between them are ignored.
// -----------------------------------------------------------------------------
void sevb_test_multiple_reactions_accumulate()
{
    std::cerr << "\n[TEST] sevb_test_multiple_reactions_accumulate\n"
              << "  Source file: src/system_setup/set_exclude_volume_bound.cpp\n"
              << "  Function:    set_exclude_volume_bound\n"
              << "  Scenario:    A(iface0) binds B(iface0) and also C(iface1);\n"
              << "               a third, opted-out, reaction sits between them.\n"
              << "  Criteria:    A.iface0 accumulates two entries in reaction\n"
              << "               order, the skipped reaction contributes nothing,\n"
              << "               and template C's uninvolved interface stays empty.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(sevb_make_mol_template(0, "A", 2));
    molTemplateList.push_back(sevb_make_mol_template(1, "B", 2));
    molTemplateList.push_back(sevb_make_mol_template(2, "C", 2));

    std::vector<ForwardRxn> forwardRxns;
    // (0) opted in : A.iface0 <-> B.iface0, radius 1.0, rxn index 0
    forwardRxns.push_back(sevb_make_rxn(true, 1.0, 0, 0, 0, 100, 1, 0, 101));
    // (1) opted OUT: A.iface0 <-> B.iface1, must be ignored completely
    forwardRxns.push_back(sevb_make_rxn(false, 9.9, 1, 0, 0, 102, 1, 1, 103));
    // (2) opted in : A.iface0 <-> C.iface1, radius 2.25, rxn index 2
    forwardRxns.push_back(sevb_make_rxn(true, 2.25, 2, 0, 0, 104, 2, 1, 105));

    std::cerr << "  Calling set_exclude_volume_bound with 3 reactions (2 opted in)...\n";
    set_exclude_volume_bound(forwardRxns, molTemplateList);

    const Interface& aIface0 = molTemplateList[0].interfaceList[0];
    std::cerr << "  Resulting interface contents:\n";
    sevb_dump_iface(molTemplateList[0], 0);
    sevb_dump_iface(molTemplateList[1], 0);
    sevb_dump_iface(molTemplateList[2], 1);

    // A.iface0 collected both opted-in partners, in reaction order.
    ASSERT_EQ(aIface0.excludeVolumeBoundList.size(), 2u)
        << "A.iface0 should have collected exactly two partners";
    EXPECT_EQ(aIface0.excludeVolumeBoundList[0], 1)
        << "first partner molecule of A.iface0 should be B (index 1)";
    EXPECT_EQ(aIface0.excludeVolumeBoundList[1], 2)
        << "second partner molecule of A.iface0 should be C (index 2)";

    ASSERT_EQ(aIface0.excludeVolumeBoundIfaceList.size(), 2u)
        << "A.iface0 should have collected exactly two partner interfaces";
    EXPECT_EQ(aIface0.excludeVolumeBoundIfaceList[0], 0)
        << "first partner interface should be B.iface0";
    EXPECT_EQ(aIface0.excludeVolumeBoundIfaceList[1], 1)
        << "second partner interface should be C.iface1";

    ASSERT_EQ(aIface0.excludeRadiusList.size(), 2u)
        << "A.iface0 should have collected exactly two radii";
    EXPECT_DOUBLE_EQ(aIface0.excludeRadiusList[0], 1.0)
        << "first radius should come from reaction 0";
    EXPECT_DOUBLE_EQ(aIface0.excludeRadiusList[1], 2.25)
        << "second radius should come from reaction 2";

    ASSERT_EQ(aIface0.excludeVolumeBoundReactList.size(), 2u)
        << "A.iface0 should have collected exactly two reaction indices";
    EXPECT_EQ(aIface0.excludeVolumeBoundReactList[0], 0)
        << "first recorded reaction index should be 0";
    EXPECT_EQ(aIface0.excludeVolumeBoundReactList[1], 2)
        << "second recorded reaction index should be 2 (the opted-out rxn 1 is skipped)";

    // The partner side of the skipped reaction (B.iface1) must remain empty.
    EXPECT_TRUE(molTemplateList[1].interfaceList[1].excludeVolumeBoundList.empty())
        << "B.iface1 only appeared in the opted-out reaction and must stay empty";

    // Partner sides of the opted-in reactions each have exactly one entry.
    ASSERT_EQ(molTemplateList[1].interfaceList[0].excludeVolumeBoundList.size(), 1u)
        << "B.iface0 should have one recorded partner";
    EXPECT_EQ(molTemplateList[1].interfaceList[0].excludeVolumeBoundList[0], 0)
        << "B.iface0 partner should be A (index 0)";
    ASSERT_EQ(molTemplateList[2].interfaceList[1].excludeVolumeBoundList.size(), 1u)
        << "C.iface1 should have one recorded partner";
    EXPECT_EQ(molTemplateList[2].interfaceList[1].excludeVolumeBoundList[0], 0)
        << "C.iface1 partner should be A (index 0)";

    // C's other interface never participated.
    EXPECT_TRUE(molTemplateList[2].interfaceList[0].excludeVolumeBoundList.empty())
        << "C.iface0 was never used and must stay empty";

    // All three templates that participated in an opted-in reaction are flagged.
    EXPECT_TRUE(molTemplateList[0].excludeVolumeBound) << "A should be flagged";
    EXPECT_TRUE(molTemplateList[1].excludeVolumeBound) << "B should be flagged";
    EXPECT_TRUE(molTemplateList[2].excludeVolumeBound) << "C should be flagged";
}

// -----------------------------------------------------------------------------
// Test 5: a homodimerization reaction (same molecule, same interface on both
//         sides) writes to the very same interface twice.
// -----------------------------------------------------------------------------
void sevb_test_self_reaction_double_entry()
{
    std::cerr << "\n[TEST] sevb_test_self_reaction_double_entry\n"
              << "  Source file: src/system_setup/set_exclude_volume_bound.cpp\n"
              << "  Function:    set_exclude_volume_bound\n"
              << "  Scenario:    homodimer A(iface0) + A(iface0), opted in,\n"
              << "               bindRadius=4.0, relRxnIndex=5.\n"
              << "  Criteria:    because both pushes target the same interface,\n"
              << "               every list on A.iface0 gains TWO identical\n"
              << "               entries and the template flag is true.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(sevb_make_mol_template(0, "A", 1));

    std::vector<ForwardRxn> forwardRxns;
    // Both reactants are molecule 0, relative interface 0.
    forwardRxns.push_back(sevb_make_rxn(true, 4.0, 5, 0, 0, 20, 0, 0, 20));

    std::cerr << "  Calling set_exclude_volume_bound on a symmetric reaction...\n";
    set_exclude_volume_bound(forwardRxns, molTemplateList);

    const Interface& aIface0 = molTemplateList[0].interfaceList[0];
    std::cerr << "  Resulting interface contents:\n";
    sevb_dump_iface(molTemplateList[0], 0);

    ASSERT_EQ(aIface0.excludeVolumeBoundList.size(), 2u)
        << "self-reaction writes the partner molecule index twice";
    EXPECT_EQ(aIface0.excludeVolumeBoundList[0], 0)
        << "first entry should be molecule 0 (itself)";
    EXPECT_EQ(aIface0.excludeVolumeBoundList[1], 0)
        << "second entry should also be molecule 0 (itself)";

    ASSERT_EQ(aIface0.excludeVolumeBoundIfaceList.size(), 2u)
        << "self-reaction writes the partner interface index twice";
    EXPECT_EQ(aIface0.excludeVolumeBoundIfaceList[0], 0)
        << "first partner interface should be 0";
    EXPECT_EQ(aIface0.excludeVolumeBoundIfaceList[1], 0)
        << "second partner interface should be 0";

    ASSERT_EQ(aIface0.excludeRadiusList.size(), 2u)
        << "self-reaction writes the bind radius twice";
    EXPECT_DOUBLE_EQ(aIface0.excludeRadiusList[0], 4.0)
        << "first radius should be 4.0";
    EXPECT_DOUBLE_EQ(aIface0.excludeRadiusList[1], 4.0)
        << "second radius should be 4.0";

    ASSERT_EQ(aIface0.excludeVolumeBoundReactList.size(), 2u)
        << "self-reaction writes the reaction index twice";
    EXPECT_EQ(aIface0.excludeVolumeBoundReactList[0], 5)
        << "first reaction index should be 5";
    EXPECT_EQ(aIface0.excludeVolumeBoundReactList[1], 5)
        << "second reaction index should be 5";

    EXPECT_TRUE(molTemplateList[0].excludeVolumeBound)
        << "MolTemplate A should be flagged excludeVolumeBound";
}

// -----------------------------------------------------------------------------
// Test 6: the function appends to, and never clears, pre-existing data, and it
//         never resets an already-true excludeVolumeBound flag.
// -----------------------------------------------------------------------------
void sevb_test_appends_to_existing_data()
{
    std::cerr << "\n[TEST] sevb_test_appends_to_existing_data\n"
              << "  Source file: src/system_setup/set_exclude_volume_bound.cpp\n"
              << "  Function:    set_exclude_volume_bound\n"
              << "  Scenario:    A.iface0 already carries a stale entry (from an\n"
              << "               earlier call) before one new opted-in reaction.\n"
              << "  Criteria:    the pre-existing entry survives at position 0 and\n"
              << "               the new data is appended at position 1.\n";

    std::vector<MolTemplate> molTemplateList;
    molTemplateList.push_back(sevb_make_mol_template(0, "A", 1));
    molTemplateList.push_back(sevb_make_mol_template(1, "B", 1));

    // Seed A.iface0 with pre-existing content, as if the routine had already
    // been run once with a different reaction list.
    molTemplateList[0].interfaceList[0].excludeVolumeBoundList.push_back(99);
    molTemplateList[0].interfaceList[0].excludeVolumeBoundIfaceList.push_back(88);
    molTemplateList[0].interfaceList[0].excludeRadiusList.push_back(0.125);
    molTemplateList[0].interfaceList[0].excludeVolumeBoundReactList.push_back(77);
    molTemplateList[0].excludeVolumeBound = true; // already flagged

    std::vector<ForwardRxn> forwardRxns;
    forwardRxns.push_back(sevb_make_rxn(true, 6.5, 3, 0, 0, 200, 1, 0, 201));

    std::cerr << "  Calling set_exclude_volume_bound on pre-populated templates...\n";
    set_exclude_volume_bound(forwardRxns, molTemplateList);

    const Interface& aIface0 = molTemplateList[0].interfaceList[0];
    std::cerr << "  Resulting interface contents:\n";
    sevb_dump_iface(molTemplateList[0], 0);

    ASSERT_EQ(aIface0.excludeVolumeBoundList.size(), 2u)
        << "the stale entry plus the new one make two entries";
    EXPECT_EQ(aIface0.excludeVolumeBoundList[0], 99)
        << "pre-existing molecule index 99 must be preserved";
    EXPECT_EQ(aIface0.excludeVolumeBoundList[1], 1)
        << "the newly appended partner molecule should be B (index 1)";

    ASSERT_EQ(aIface0.excludeVolumeBoundIfaceList.size(), 2u)
        << "the stale interface entry plus the new one make two entries";
    EXPECT_EQ(aIface0.excludeVolumeBoundIfaceList[0], 88)
        << "pre-existing interface index 88 must be preserved";
    EXPECT_EQ(aIface0.excludeVolumeBoundIfaceList[1], 0)
        << "the newly appended partner interface should be 0";

    ASSERT_EQ(aIface0.excludeRadiusList.size(), 2u)
        << "the stale radius plus the new one make two entries";
    EXPECT_DOUBLE_EQ(aIface0.excludeRadiusList[0], 0.125)
        << "pre-existing radius 0.125 must be preserved";
    EXPECT_DOUBLE_EQ(aIface0.excludeRadiusList[1], 6.5)
        << "the newly appended radius should be 6.5";

    ASSERT_EQ(aIface0.excludeVolumeBoundReactList.size(), 2u)
        << "the stale reaction index plus the new one make two entries";
    EXPECT_EQ(aIface0.excludeVolumeBoundReactList[0], 77)
        << "pre-existing reaction index 77 must be preserved";
    EXPECT_EQ(aIface0.excludeVolumeBoundReactList[1], 3)
        << "the newly appended reaction index should be 3";

    EXPECT_TRUE(molTemplateList[0].excludeVolumeBound)
        << "an already-true flag must remain true";
    EXPECT_TRUE(molTemplateList[1].excludeVolumeBound)
        << "the partner template B must now be flagged as well";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named helper runs inside its own TEST so that a
// failure in one scenario does not prevent the others from executing.
// -----------------------------------------------------------------------------
TEST(SetExcludeVolumeBound, EmptyReactionList) { sevb_test_empty_reaction_list(); }
TEST(SetExcludeVolumeBound, FlagFalseIsSkipped) { sevb_test_flag_false_is_skipped(); }
TEST(SetExcludeVolumeBound, SingleReactionCrossRegisters) { sevb_test_single_reaction_cross_registers(); }
TEST(SetExcludeVolumeBound, MultipleReactionsAccumulate) { sevb_test_multiple_reactions_accumulate(); }
TEST(SetExcludeVolumeBound, SelfReactionDoubleEntry) { sevb_test_self_reaction_double_entry(); }
TEST(SetExcludeVolumeBound, AppendsToExistingData) { sevb_test_appends_to_existing_data(); }