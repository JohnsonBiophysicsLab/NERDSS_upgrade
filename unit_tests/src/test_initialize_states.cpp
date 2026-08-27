/*! \file test_initialize_states.cpp
 *
 * ### Unit test for src/system_setup/initialize_states.cpp
 *
 * Function under test:
 *
 *     void initialize_states(std::vector<Molecule>& moleculeList,
 *                            std::vector<MolTemplate>& molTemplateList,
 *                            Membrane& membraneObject)
 *
 * What the function does (from reading the implementation):
 *
 *  - For every MolTemplate it looks at `startingNumState`, a small container
 *    holding `totalCopyNumbers`, `numberEachState` and `nameEachState`.
 *    `nameEachState` entries look like "a~U" or, for multi-interface
 *    definitions, "b1~P,b2~Y".  Each entry is parsed character by character:
 *      * alphanumeric characters are appended to a buffer (the interface name),
 *      * '~' terminates the interface name, and the *next* character is the
 *        state identifier, converted to UPPER CASE,
 *      * ',' clears the buffer so a second interface can be read,
 *      * anything else is a fatal parse error (the code calls exit(1), so this
 *        test deliberately never feeds it an invalid character).
 *
 *  - If the template is an implicit lipid, the parsed copy numbers are written
 *    into `membraneObject.numberOfFreeLipidsEachState[stateIndex]`, where
 *    stateIndex is the position of the matching state in
 *    `interfaceList[0].stateList`.  The moleculeList is *not* touched.
 *    NOTE: `numberOfFreeLipidsEachState` is indexed without bounds checks, so
 *    the test always pre-sizes it.
 *
 *  - Otherwise the first `totalCopyNumbers` molecules of that molecule type
 *    (starting at the first molecule in moleculeList whose molTypeIndex
 *    matches) have their interface state fields (index, relIndex, stateIden,
 *    stateIndex) overwritten.  Molecules are assigned to groups by a running
 *    cumulative sum over `numberEachState`.
 *
 * Preconditions the implementation silently assumes (and which this test
 * therefore always satisfies, rather than triggering undefined behaviour):
 *   - at least one molecule of the templated type exists in moleculeList,
 *   - sum(numberEachState) >= totalCopyNumbers,
 *   - every parsed interface name / state identifier actually exists on the
 *     template,
 *   - membraneObject.numberOfFreeLipidsEachState is large enough.
 */

#include "system_setup/system_setup.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <string>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
// Local helpers (anonymous namespace + unique prefix so nothing collides with
// other translation units in the combined gtest binary).
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a fully initialized MolTemplate Interface.
 *
 * \param[in] name     interface name, e.g. "a"
 * \param[in] relIndex the interface's relative index inside the template
 * \param[in] states   list of (state identifier, absolute state index) pairs
 */
Interface initstates_make_interface(const std::string& name, int relIndex,
    const std::vector<std::pair<char, int>>& states)
{
    Interface iface {};
    iface.name = name;
    iface.index = relIndex;
    iface.iCoord = Coord { 0.0, 0.0, 0.0 };
    for (const auto& oneState : states) {
        // Interface::State(char iden, int index)
        Interface::State state(oneState.first, oneState.second);
        state.ifaceAndStateName = name + "~" + oneState.first;
        iface.stateList.push_back(state);
    }
    return iface;
}

/*! \brief Build a MolTemplate with everything initialize_states() reads. */
MolTemplate initstates_make_template(const std::string& molName, int molTypeIndex,
    const std::vector<Interface>& ifaces, bool isImplicitLipid = false)
{
    MolTemplate tmpl {};
    tmpl.molName = molName;
    tmpl.molTypeIndex = molTypeIndex;
    tmpl.interfaceList = ifaces;
    tmpl.isImplicitLipid = isImplicitLipid;
    tmpl.isLipid = isImplicitLipid;
    tmpl.mass = 1.0;
    tmpl.radius = 1.0;
    tmpl.copies = 0;
    tmpl.D = Coord { 1.0, 1.0, 1.0 };
    tmpl.Dr = Coord { 0.01, 0.01, 0.01 };
    // startingNumState is intentionally left empty; each test fills it in.
    return tmpl;
}

/*! \brief Build a Molecule whose interfaces sit in the *first* (default) state.
 *
 * This mirrors how Molecule::create_random_coords() initializes a real
 * molecule: every interface gets stateList[0]'s identity and absolute index.
 *
 * \param[in] index the molecule's index inside moleculeList (must equal the
 *                  position in the vector, because initialize_states() does
 *                  pointer arithmetic on it)
 * \param[in] tmpl  the template this molecule is an instance of
 */
Molecule initstates_make_molecule(int index, const MolTemplate& tmpl)
{
    Molecule mol {};
    mol.index = index;
    mol.id = index;
    mol.myComIndex = index;
    mol.molTypeIndex = tmpl.molTypeIndex;
    mol.mass = tmpl.mass;
    mol.isLipid = tmpl.isLipid;
    mol.isImplicitLipid = tmpl.isImplicitLipid;
    mol.comCoord = Coord { 0.0, 0.0, 0.0 };

    for (const auto& tIface : tmpl.interfaceList) {
        Molecule::Iface iface {};
        iface.coord = Coord { 0.0, 0.0, 0.0 };
        iface.relIndex = tIface.index;
        iface.molTypeIndex = tmpl.molTypeIndex;
        iface.stateIndex = 0; // default: first state in the list
        iface.index = tIface.stateList[0].index;
        iface.stateIden = tIface.stateList[0].iden;
        iface.isBound = false;
        mol.interfaceList.push_back(iface);
    }
    return mol;
}

/*! \brief Convenience: print one molecule's interface state to stderr. */
void initstates_dump_molecule(const Molecule& mol)
{
    std::cerr << "      mol " << mol.index << " :";
    for (const auto& iface : mol.interfaceList) {
        std::cerr << " [rel=" << iface.relIndex << " abs=" << iface.index
                  << " iden=" << iface.stateIden << " stateIdx=" << iface.stateIndex << "]";
    }
    std::cerr << '\n';
}

} // anonymous namespace

// -----------------------------------------------------------------------------
// Test 1: a template with an empty startingNumState must be a complete no-op.
// -----------------------------------------------------------------------------
void test_initstates_no_starting_states()
{
    std::cerr << "\n[TEST] test_initstates_no_starting_states\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states()\n"
              << "  Scenario:    MolTemplate has states available but the parsed\n"
              << "               startingNumState block is empty (nameEachState is empty).\n"
              << "  Pass:        molecules keep their default state and the membrane\n"
              << "               free-lipid counters are untouched.\n";

    // Template "A" with one interface "a" that can be U (abs 0) or P (abs 1).
    MolTemplate tmplA = initstates_make_template(
        "A", 0, { initstates_make_interface("a", 0, { { 'U', 0 }, { 'P', 1 } }) });
    // startingNumState left empty on purpose -> tmpLength == 0 -> nothing happens.

    std::vector<MolTemplate> molTemplateList { tmplA };

    std::vector<Molecule> moleculeList;
    for (int i = 0; i < 3; ++i)
        moleculeList.push_back(initstates_make_molecule(i, tmplA));

    Membrane membraneObject {};
    membraneObject.nStates = 2;
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states() with 3 molecules, no starting states...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    for (const auto& mol : moleculeList) {
        initstates_dump_molecule(mol);
        EXPECT_EQ(mol.interfaceList[0].stateIden, 'U')
            << "molecule " << mol.index << " should still be in the default state 'U'";
        EXPECT_EQ(mol.interfaceList[0].stateIndex, 0)
            << "molecule " << mol.index << " default relative state index should still be 0";
        EXPECT_EQ(mol.interfaceList[0].index, 0)
            << "molecule " << mol.index << " default absolute state index should still be 0";
    }

    // The membrane counters belong to the implicit-lipid path only.
    ASSERT_EQ(membraneObject.numberOfFreeLipidsEachState.size(), 2u)
        << "membrane counter vector should not be resized";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 0)
        << "membrane free lipid counter 0 should be untouched";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[1], 0)
        << "membrane free lipid counter 1 should be untouched";
}

// -----------------------------------------------------------------------------
// Test 2: a single interface split across two states (1 x U, 3 x P).
// -----------------------------------------------------------------------------
void test_initstates_single_interface_split()
{
    std::cerr << "\n[TEST] test_initstates_single_interface_split\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (explicit / non-implicit branch)\n"
              << "  Scenario:    nameEachState = {\"a~U\", \"a~P\"},\n"
              << "               numberEachState = {1, 3}, totalCopyNumbers = 4.\n"
              << "  Pass:        molecule 0 keeps state 'U'; molecules 1,2,3 become 'P'\n"
              << "               with absolute index 1 and relative state index 1.\n";

    MolTemplate tmplA = initstates_make_template(
        "A", 0, { initstates_make_interface("a", 0, { { 'U', 0 }, { 'P', 1 } }) });
    tmplA.startingNumState.totalCopyNumbers = 4;
    tmplA.startingNumState.nameEachState = { "a~U", "a~P" };
    tmplA.startingNumState.numberEachState = { 1, 3 };

    std::vector<MolTemplate> molTemplateList { tmplA };

    std::vector<Molecule> moleculeList;
    for (int i = 0; i < 4; ++i)
        moleculeList.push_back(initstates_make_molecule(i, tmplA));

    Membrane membraneObject {};
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    std::cerr << "  Resulting molecule states:\n";
    for (const auto& mol : moleculeList)
        initstates_dump_molecule(mol);

    // Molecule 0 falls into the first group ("a~U").
    EXPECT_EQ(moleculeList[0].interfaceList[0].stateIden, 'U')
        << "the first molecule should be assigned state 'U'";
    EXPECT_EQ(moleculeList[0].interfaceList[0].stateIndex, 0)
        << "'U' is stateList[0], so the relative state index should be 0";
    EXPECT_EQ(moleculeList[0].interfaceList[0].index, 0)
        << "'U' has absolute state index 0";
    EXPECT_EQ(moleculeList[0].interfaceList[0].relIndex, 0)
        << "relIndex must be set to the interface's relative index (0)";

    // Molecules 1..3 fall into the second group ("a~P").
    for (int i = 1; i < 4; ++i) {
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIden, 'P')
            << "molecule " << i << " should be assigned state 'P'";
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIndex, 1)
            << "'P' is stateList[1], so the relative state index should be 1";
        EXPECT_EQ(moleculeList[i].interfaceList[0].index, 1)
            << "'P' has absolute state index 1";
        EXPECT_EQ(moleculeList[i].interfaceList[0].relIndex, 0)
            << "relIndex must be set to the interface's relative index (0)";
    }
}

// -----------------------------------------------------------------------------
// Test 3: a group with zero copies must simply be skipped by the cumulative
//         counting loop.
// -----------------------------------------------------------------------------
void test_initstates_zero_count_group()
{
    std::cerr << "\n[TEST] test_initstates_zero_count_group\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (cumulative group selection)\n"
              << "  Scenario:    numberEachState = {0, 3} -> the first group holds no\n"
              << "               molecules at all.\n"
              << "  Pass:        every one of the 3 molecules lands in the second\n"
              << "               group and becomes 'P'.\n";

    MolTemplate tmplA = initstates_make_template(
        "A", 0, { initstates_make_interface("a", 0, { { 'U', 0 }, { 'P', 1 } }) });
    tmplA.startingNumState.totalCopyNumbers = 3;
    tmplA.startingNumState.nameEachState = { "a~U", "a~P" };
    tmplA.startingNumState.numberEachState = { 0, 3 };

    std::vector<MolTemplate> molTemplateList { tmplA };

    std::vector<Molecule> moleculeList;
    for (int i = 0; i < 3; ++i)
        moleculeList.push_back(initstates_make_molecule(i, tmplA));

    Membrane membraneObject {};
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    for (const auto& mol : moleculeList) {
        initstates_dump_molecule(mol);
        EXPECT_EQ(mol.interfaceList[0].stateIden, 'P')
            << "molecule " << mol.index << " should skip the empty group and become 'P'";
        EXPECT_EQ(mol.interfaceList[0].stateIndex, 1)
            << "molecule " << mol.index << " relative state index should be 1";
        EXPECT_EQ(mol.interfaceList[0].index, 1)
            << "molecule " << mol.index << " absolute state index should be 1";
    }
}

// -----------------------------------------------------------------------------
// Test 4: a comma-separated entry sets the state of two interfaces at once.
// -----------------------------------------------------------------------------
void test_initstates_multiple_interfaces()
{
    std::cerr << "\n[TEST] test_initstates_multiple_interfaces\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (comma parsing)\n"
              << "  Scenario:    nameEachState = {\"b1~P,b2~Y\"}, numberEachState = {2}.\n"
              << "  Pass:        both molecules get interface 0 -> 'P' (abs 1) and\n"
              << "               interface 1 -> 'Y' (abs 3).\n";

    MolTemplate tmplB = initstates_make_template("B", 0,
        { initstates_make_interface("b1", 0, { { 'U', 0 }, { 'P', 1 } }),
            initstates_make_interface("b2", 1, { { 'X', 2 }, { 'Y', 3 } }) });
    tmplB.startingNumState.totalCopyNumbers = 2;
    tmplB.startingNumState.nameEachState = { "b1~P,b2~Y" };
    tmplB.startingNumState.numberEachState = { 2 };

    std::vector<MolTemplate> molTemplateList { tmplB };

    std::vector<Molecule> moleculeList;
    for (int i = 0; i < 2; ++i)
        moleculeList.push_back(initstates_make_molecule(i, tmplB));

    Membrane membraneObject {};
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    for (const auto& mol : moleculeList) {
        initstates_dump_molecule(mol);

        // Interface "b1" -> state 'P'
        EXPECT_EQ(mol.interfaceList[0].stateIden, 'P')
            << "molecule " << mol.index << " interface b1 should be 'P'";
        EXPECT_EQ(mol.interfaceList[0].stateIndex, 1)
            << "molecule " << mol.index << " interface b1 relative state index should be 1";
        EXPECT_EQ(mol.interfaceList[0].index, 1)
            << "molecule " << mol.index << " interface b1 absolute state index should be 1";
        EXPECT_EQ(mol.interfaceList[0].relIndex, 0)
            << "molecule " << mol.index << " interface b1 relIndex should be 0";

        // Interface "b2" -> state 'Y'
        EXPECT_EQ(mol.interfaceList[1].stateIden, 'Y')
            << "molecule " << mol.index << " interface b2 should be 'Y'";
        EXPECT_EQ(mol.interfaceList[1].stateIndex, 1)
            << "molecule " << mol.index << " interface b2 relative state index should be 1";
        EXPECT_EQ(mol.interfaceList[1].index, 3)
            << "molecule " << mol.index << " interface b2 absolute state index should be 3";
        EXPECT_EQ(mol.interfaceList[1].relIndex, 1)
            << "molecule " << mol.index << " interface b2 relIndex should be 1";
    }
}

// -----------------------------------------------------------------------------
// Test 5: the state character is upper-cased during parsing, so "a~p" must
//         match the template state 'P'.
// -----------------------------------------------------------------------------
void test_initstates_lowercase_state_is_uppercased()
{
    std::cerr << "\n[TEST] test_initstates_lowercase_state_is_uppercased\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (std::toupper on the state char)\n"
              << "  Scenario:    nameEachState = {\"a~p\"} (lower case state letter).\n"
              << "  Pass:        molecules are matched against the template state 'P'.\n";

    MolTemplate tmplA = initstates_make_template(
        "A", 0, { initstates_make_interface("a", 0, { { 'U', 0 }, { 'P', 1 } }) });
    tmplA.startingNumState.totalCopyNumbers = 2;
    tmplA.startingNumState.nameEachState = { "a~p" };
    tmplA.startingNumState.numberEachState = { 2 };

    std::vector<MolTemplate> molTemplateList { tmplA };

    std::vector<Molecule> moleculeList;
    for (int i = 0; i < 2; ++i)
        moleculeList.push_back(initstates_make_molecule(i, tmplA));

    Membrane membraneObject {};
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    for (const auto& mol : moleculeList) {
        initstates_dump_molecule(mol);
        EXPECT_EQ(mol.interfaceList[0].stateIden, 'P')
            << "lower-case 'p' in the input should resolve to template state 'P'";
        EXPECT_EQ(mol.interfaceList[0].stateIndex, 1)
            << "the resolved state should be stateList[1]";
        EXPECT_EQ(mol.interfaceList[0].index, 1)
            << "the resolved state has absolute index 1";
    }
}

// -----------------------------------------------------------------------------
// Test 6: molecules of a different type that appear *before* the templated
//         molecules in moleculeList must be left alone (firstIndex offset).
// -----------------------------------------------------------------------------
void test_initstates_offset_molecule_indices()
{
    std::cerr << "\n[TEST] test_initstates_offset_molecule_indices\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (firstIndex search + offset)\n"
              << "  Scenario:    moleculeList = [A0, A1, B2, B3, B4]; only template B\n"
              << "               declares starting states.\n"
              << "  Pass:        molecules 0 and 1 keep their default 'U' state, while\n"
              << "               molecules 2,3,4 become 'P'.\n";

    MolTemplate tmplA = initstates_make_template(
        "A", 0, { initstates_make_interface("a", 0, { { 'U', 0 }, { 'P', 1 } }) });
    // Template A intentionally has no startingNumState block.

    MolTemplate tmplB = initstates_make_template(
        "B", 1, { initstates_make_interface("b", 0, { { 'U', 2 }, { 'P', 3 } }) });
    tmplB.startingNumState.totalCopyNumbers = 3;
    tmplB.startingNumState.nameEachState = { "b~P" };
    tmplB.startingNumState.numberEachState = { 3 };

    std::vector<MolTemplate> molTemplateList { tmplA, tmplB };

    std::vector<Molecule> moleculeList;
    moleculeList.push_back(initstates_make_molecule(0, tmplA));
    moleculeList.push_back(initstates_make_molecule(1, tmplA));
    moleculeList.push_back(initstates_make_molecule(2, tmplB));
    moleculeList.push_back(initstates_make_molecule(3, tmplB));
    moleculeList.push_back(initstates_make_molecule(4, tmplB));

    Membrane membraneObject {};
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    std::cerr << "  Resulting molecule states:\n";
    for (const auto& mol : moleculeList)
        initstates_dump_molecule(mol);

    // Type A molecules are untouched.
    for (int i = 0; i < 2; ++i) {
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIden, 'U')
            << "type-A molecule " << i << " must not be modified";
        EXPECT_EQ(moleculeList[i].interfaceList[0].index, 0)
            << "type-A molecule " << i << " keeps absolute state index 0";
    }

    // Type B molecules all become 'P' (absolute index 3).
    for (int i = 2; i < 5; ++i) {
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIden, 'P')
            << "type-B molecule " << i << " should be set to 'P'";
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIndex, 1)
            << "type-B molecule " << i << " relative state index should be 1";
        EXPECT_EQ(moleculeList[i].interfaceList[0].index, 3)
            << "type-B molecule " << i << " absolute state index should be 3";
    }
}

// -----------------------------------------------------------------------------
// Test 7: only the first `totalCopyNumbers` molecules of a type are touched.
// -----------------------------------------------------------------------------
void test_initstates_partial_total_copy_numbers()
{
    std::cerr << "\n[TEST] test_initstates_partial_total_copy_numbers\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (loop bound = totalCopyNumbers)\n"
              << "  Scenario:    4 molecules of type A exist but totalCopyNumbers = 2.\n"
              << "  Pass:        molecules 0,1 become 'P'; molecules 2,3 stay 'U'.\n";

    MolTemplate tmplA = initstates_make_template(
        "A", 0, { initstates_make_interface("a", 0, { { 'U', 0 }, { 'P', 1 } }) });
    tmplA.startingNumState.totalCopyNumbers = 2;
    tmplA.startingNumState.nameEachState = { "a~P" };
    tmplA.startingNumState.numberEachState = { 2 };

    std::vector<MolTemplate> molTemplateList { tmplA };

    std::vector<Molecule> moleculeList;
    for (int i = 0; i < 4; ++i)
        moleculeList.push_back(initstates_make_molecule(i, tmplA));

    Membrane membraneObject {};
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    std::cerr << "  Resulting molecule states:\n";
    for (const auto& mol : moleculeList)
        initstates_dump_molecule(mol);

    for (int i = 0; i < 2; ++i) {
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIden, 'P')
            << "molecule " << i << " is inside totalCopyNumbers and should be 'P'";
        EXPECT_EQ(moleculeList[i].interfaceList[0].index, 1)
            << "molecule " << i << " absolute state index should be 1";
    }
    for (int i = 2; i < 4; ++i) {
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIden, 'U')
            << "molecule " << i << " is beyond totalCopyNumbers and must keep 'U'";
        EXPECT_EQ(moleculeList[i].interfaceList[0].index, 0)
            << "molecule " << i << " absolute state index should still be 0";
    }
}

// -----------------------------------------------------------------------------
// Test 8: implicit-lipid templates write into the membrane free-lipid counters
//         and leave moleculeList completely alone.
// -----------------------------------------------------------------------------
void test_initstates_implicit_lipid_membrane_counts()
{
    std::cerr << "\n[TEST] test_initstates_implicit_lipid_membrane_counts\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (implicit-lipid branch)\n"
              << "  Scenario:    implicit lipid with states A/B/C, starting numbers\n"
              << "               {\"il~B\": 7, \"il~C\": 5}.\n"
              << "  Pass:        numberOfFreeLipidsEachState becomes {0, 7, 5} and the\n"
              << "               implicit-lipid molecule itself is untouched.\n";

    MolTemplate tmplIL = initstates_make_template("il", 0,
        { initstates_make_interface("il", 0, { { 'A', 0 }, { 'B', 1 }, { 'C', 2 } }) },
        /*isImplicitLipid=*/true);
    tmplIL.startingNumState.totalCopyNumbers = 12;
    tmplIL.startingNumState.nameEachState = { "il~B", "il~C" };
    tmplIL.startingNumState.numberEachState = { 7, 5 };

    std::vector<MolTemplate> molTemplateList { tmplIL };

    // A single implicit-lipid "molecule" placeholder, as the real code creates.
    std::vector<Molecule> moleculeList { initstates_make_molecule(0, tmplIL) };

    Membrane membraneObject {};
    membraneObject.implicitLipid = true;
    membraneObject.nStates = 3;
    // MUST be pre-sized: the implementation indexes it without bounds checking.
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    std::cerr << "  numberOfFreeLipidsEachState = {";
    for (size_t i = 0; i < membraneObject.numberOfFreeLipidsEachState.size(); ++i)
        std::cerr << (i ? ", " : "") << membraneObject.numberOfFreeLipidsEachState[i];
    std::cerr << "}\n";

    ASSERT_EQ(membraneObject.numberOfFreeLipidsEachState.size(), 3u)
        << "the counter vector should keep its size";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 0)
        << "state 'A' was never listed, so its counter must stay 0";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[1], 7)
        << "state 'B' is stateList[1] and was given 7 copies";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[2], 5)
        << "state 'C' is stateList[2] and was given 5 copies";

    // The implicit branch never touches moleculeList.
    initstates_dump_molecule(moleculeList[0]);
    EXPECT_EQ(moleculeList[0].interfaceList[0].stateIden, 'A')
        << "the implicit-lipid molecule should keep its default state";
    EXPECT_EQ(moleculeList[0].interfaceList[0].stateIndex, 0)
        << "the implicit-lipid molecule relative state index should stay 0";
    EXPECT_EQ(moleculeList[0].interfaceList[0].index, 0)
        << "the implicit-lipid molecule absolute state index should stay 0";
}

// -----------------------------------------------------------------------------
// Test 9: several templates handled in one call - an explicit one and an
//         implicit one - each taking its own code path.
// -----------------------------------------------------------------------------
void test_initstates_multiple_templates_one_call()
{
    std::cerr << "\n[TEST] test_initstates_multiple_templates_one_call\n"
              << "  Source file: src/system_setup/initialize_states.cpp\n"
              << "  Function:    initialize_states() (loop over molTemplateList)\n"
              << "  Scenario:    templateList = [implicit lipid, explicit protein];\n"
              << "               both declare starting states.\n"
              << "  Pass:        the membrane counters AND the protein molecule states\n"
              << "               are both updated in the same call.\n";

    // Template 0: implicit lipid with two states.
    MolTemplate tmplIL = initstates_make_template("il", 0,
        { initstates_make_interface("il", 0, { { 'A', 0 }, { 'B', 1 } }) },
        /*isImplicitLipid=*/true);
    tmplIL.startingNumState.totalCopyNumbers = 10;
    tmplIL.startingNumState.nameEachState = { "il~A", "il~B" };
    tmplIL.startingNumState.numberEachState = { 4, 6 };

    // Template 1: an ordinary protein with two interfaces.
    MolTemplate tmplP = initstates_make_template("prot", 1,
        { initstates_make_interface("p1", 0, { { 'U', 2 }, { 'P', 3 } }),
            initstates_make_interface("p2", 1, { { 'X', 4 }, { 'Y', 5 } }) });
    tmplP.startingNumState.totalCopyNumbers = 3;
    tmplP.startingNumState.nameEachState = { "p1~P,p2~Y", "p1~U,p2~X" };
    tmplP.startingNumState.numberEachState = { 1, 2 };

    std::vector<MolTemplate> molTemplateList { tmplIL, tmplP };

    // moleculeList: one implicit lipid placeholder followed by three proteins.
    std::vector<Molecule> moleculeList;
    moleculeList.push_back(initstates_make_molecule(0, tmplIL));
    for (int i = 1; i < 4; ++i)
        moleculeList.push_back(initstates_make_molecule(i, tmplP));

    Membrane membraneObject {};
    membraneObject.implicitLipid = true;
    membraneObject.nStates = 2;
    membraneObject.numberOfFreeLipidsEachState = std::vector<int> { 0, 0 };

    std::cerr << "  Calling initialize_states()...\n";
    initialize_states(moleculeList, molTemplateList, membraneObject);

    // --- implicit lipid side ---
    std::cerr << "  numberOfFreeLipidsEachState = {"
              << membraneObject.numberOfFreeLipidsEachState[0] << ", "
              << membraneObject.numberOfFreeLipidsEachState[1] << "}\n";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[0], 4)
        << "implicit state 'A' should hold 4 free lipids";
    EXPECT_EQ(membraneObject.numberOfFreeLipidsEachState[1], 6)
        << "implicit state 'B' should hold 6 free lipids";

    // --- explicit protein side ---
    std::cerr << "  Resulting protein states:\n";
    for (size_t i = 1; i < moleculeList.size(); ++i)
        initstates_dump_molecule(moleculeList[i]);

    // Protein at moleculeList index 1 is the single "p1~P,p2~Y" copy.
    EXPECT_EQ(moleculeList[1].interfaceList[0].stateIden, 'P')
        << "first protein interface p1 should be 'P'";
    EXPECT_EQ(moleculeList[1].interfaceList[0].index, 3)
        << "first protein interface p1 absolute state index should be 3";
    EXPECT_EQ(moleculeList[1].interfaceList[1].stateIden, 'Y')
        << "first protein interface p2 should be 'Y'";
    EXPECT_EQ(moleculeList[1].interfaceList[1].index, 5)
        << "first protein interface p2 absolute state index should be 5";

    // Proteins at moleculeList indices 2 and 3 are the two "p1~U,p2~X" copies.
    for (int i = 2; i < 4; ++i) {
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIden, 'U')
            << "protein " << i << " interface p1 should be 'U'";
        EXPECT_EQ(moleculeList[i].interfaceList[0].stateIndex, 0)
            << "protein " << i << " interface p1 relative state index should be 0";
        EXPECT_EQ(moleculeList[i].interfaceList[0].index, 2)
            << "protein " << i << " interface p1 absolute state index should be 2";
        EXPECT_EQ(moleculeList[i].interfaceList[1].stateIden, 'X')
            << "protein " << i << " interface p2 should be 'X'";
        EXPECT_EQ(moleculeList[i].interfaceList[1].stateIndex, 0)
            << "protein " << i << " interface p2 relative state index should be 0";
        EXPECT_EQ(moleculeList[i].interfaceList[1].index, 4)
            << "protein " << i << " interface p2 absolute state index should be 4";
    }

    // The implicit-lipid placeholder molecule is still untouched.
    EXPECT_EQ(moleculeList[0].interfaceList[0].stateIden, 'A')
        << "the implicit-lipid molecule must not be modified by the explicit branch";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers: each named test_* function is run inside its own TEST so
// the framework reports individual results while still executing all of them.
// -----------------------------------------------------------------------------
TEST(InitializeStatesTest, NoStartingStatesIsNoOp) { test_initstates_no_starting_states(); }
TEST(InitializeStatesTest, SingleInterfaceSplit) { test_initstates_single_interface_split(); }
TEST(InitializeStatesTest, ZeroCountGroupSkipped) { test_initstates_zero_count_group(); }
TEST(InitializeStatesTest, MultipleInterfacesCommaSeparated) { test_initstates_multiple_interfaces(); }
TEST(InitializeStatesTest, LowercaseStateIsUppercased) { test_initstates_lowercase_state_is_uppercased(); }
TEST(InitializeStatesTest, OffsetMoleculeIndices) { test_initstates_offset_molecule_indices(); }
TEST(InitializeStatesTest, PartialTotalCopyNumbers) { test_initstates_partial_total_copy_numbers(); }
TEST(InitializeStatesTest, ImplicitLipidMembraneCounts) { test_initstates_implicit_lipid_membrane_counts(); }
TEST(InitializeStatesTest, MultipleTemplatesOneCall) { test_initstates_multiple_templates_one_call(); }