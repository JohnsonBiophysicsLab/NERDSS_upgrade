// class_MolTemplate_test.cpp
// Unit tests for classes/class_MolTemplate.hpp implementation (class_moltemplate_functions.cpp)
//
// This test exercises the Interface, Interface::State, and MolTemplate classes,
// including their constructors and member functions. Verbose console output is
// provided so that when the test suite runs the reader can follow which function
// in which source file is being tested and what each test is verifying.

#include "classes/class_MolTemplate.hpp"
#include "parser/parser_functions.hpp"

#include <gtest/gtest.h>
#include <iostream>
#include <string>
#include <vector>

// ---------------------------------------------------------------------------
// Helper: print a small banner describing what is being tested.
// ---------------------------------------------------------------------------
static void print_test_banner(const std::string& funcName, const std::string& fileName)
{
    std::cerr << "\n=============================================================\n";
    std::cerr << "[TEST] Testing function: " << funcName << "\n";
    std::cerr << "[TEST] Source file     : " << fileName << "\n";
    std::cerr << "=============================================================\n";
}

// ---------------------------------------------------------------------------
// Test Interface::State constructors.
// There are several overloads: (index), (name, iden, index), (name, index).
// We verify that the members are set correctly for each variant.
// ---------------------------------------------------------------------------
void test_State_constructors()
{
    print_test_banner("Interface::State constructors", "class_moltemplate_functions.cpp");

    // Constructor taking only an index.
    std::cerr << "[INFO] Creating State with index-only constructor (index = 5).\n";
    Interface::State s1(5);
    // Expect the index to be stored.
    EXPECT_EQ(s1.index, 5) << "State(int index) should set index to 5";
    std::cerr << "[INFO] s1.index = " << s1.index << "\n";

    // Constructor taking name, identity char, and index.
    std::cerr << "[INFO] Creating State with (name, iden, index) constructor.\n";
    Interface::State s2("site~A", 'A', 2);
    EXPECT_EQ(s2.ifaceAndStateName, "site~A") << "ifaceAndStateName should be set";
    EXPECT_EQ(s2.iden, 'A') << "iden should be 'A'";
    EXPECT_EQ(s2.index, 2) << "index should be 2";
    std::cerr << "[INFO] s2.ifaceAndStateName = " << s2.ifaceAndStateName
              << ", s2.iden = " << s2.iden << ", s2.index = " << s2.index << "\n";

    // Constructor taking name and index only.
    std::cerr << "[INFO] Creating State with (name, index) constructor.\n";
    Interface::State s3("site", 7);
    EXPECT_EQ(s3.ifaceAndStateName, "site") << "ifaceAndStateName should be 'site'";
    EXPECT_EQ(s3.index, 7) << "index should be 7";
    std::cerr << "[INFO] s3.ifaceAndStateName = " << s3.ifaceAndStateName
              << ", s3.index = " << s3.index << "\n";
}

// ---------------------------------------------------------------------------
// Test Interface constructors and set_ifaceAndStateNames().
// The (name, states, coord) constructor calls set_ifaceAndStateNames() which
// builds the combined "name~iden" string for states when there are more than one.
// ---------------------------------------------------------------------------
void test_Interface_constructors_and_naming()
{
    print_test_banner("Interface constructors + set_ifaceAndStateNames", "class_moltemplate_functions.cpp");

    // Simple (name, coord) constructor.
    std::cerr << "[INFO] Creating Interface with (name, Coord) constructor.\n";
    Coord c(1.0, 2.0, 3.0);
    Interface iface1("bindingSite", c);
    EXPECT_EQ(iface1.name, "bindingSite") << "Interface name should be set";
    EXPECT_DOUBLE_EQ(iface1.iCoord.x, 1.0) << "iCoord.x should be 1.0";
    EXPECT_DOUBLE_EQ(iface1.iCoord.y, 2.0) << "iCoord.y should be 2.0";
    EXPECT_DOUBLE_EQ(iface1.iCoord.z, 3.0) << "iCoord.z should be 3.0";
    std::cerr << "[INFO] iface1.name = " << iface1.name << "\n";

    // (name, states, coord) constructor with multiple states triggers name generation.
    std::cerr << "[INFO] Creating Interface with multiple states to test set_ifaceAndStateNames().\n";
    std::vector<Interface::State> states;
    states.emplace_back("", 'A', 0); // state A
    states.emplace_back("", 'B', 1); // state B
    Interface iface2("phos", states, c);

    // With more than one state and a valid iden, each state's ifaceAndStateName
    // should be "name~iden".
    ASSERT_EQ(iface2.stateList.size(), 2u) << "iface2 should contain two states";
    EXPECT_EQ(iface2.stateList[0].ifaceAndStateName, "phos~A")
        << "State 0 name should be combined as 'phos~A'";
    EXPECT_EQ(iface2.stateList[1].ifaceAndStateName, "phos~B")
        << "State 1 name should be combined as 'phos~B'";
    std::cerr << "[INFO] iface2.stateList[0].ifaceAndStateName = "
              << iface2.stateList[0].ifaceAndStateName << "\n";
    std::cerr << "[INFO] iface2.stateList[1].ifaceAndStateName = "
              << iface2.stateList[1].ifaceAndStateName << "\n";

    // With a single state the name should NOT be modified.
    std::cerr << "[INFO] Creating Interface with a single state (no name generation expected).\n";
    std::vector<Interface::State> singleState;
    singleState.emplace_back("origName", 'X', 0);
    Interface iface3("single", singleState, c);
    EXPECT_EQ(iface3.stateList[0].ifaceAndStateName, "origName")
        << "Single-state interface should not modify ifaceAndStateName";
    std::cerr << "[INFO] iface3.stateList[0].ifaceAndStateName = "
              << iface3.stateList[0].ifaceAndStateName << "\n";
}

// ---------------------------------------------------------------------------
// Test MolTemplate constructors.
// Two constructors exist: (molName, interfaceList) and (comCoord, interfaceList).
// ---------------------------------------------------------------------------
void test_MolTemplate_constructors()
{
    print_test_banner("MolTemplate constructors", "class_moltemplate_functions.cpp");

    // Build a small interface list.
    Coord c(0.0, 0.0, 0.0);
    std::vector<Interface> ifaceList;
    ifaceList.emplace_back("iface0", c);

    // Constructor with molecule name.
    std::cerr << "[INFO] Creating MolTemplate with (molName, interfaceList) constructor.\n";
    MolTemplate mt1("myMolecule", ifaceList);
    EXPECT_EQ(mt1.molName, "myMolecule") << "molName should be 'myMolecule'";
    EXPECT_EQ(mt1.interfaceList.size(), 1u) << "interfaceList should have 1 element";
    std::cerr << "[INFO] mt1.molName = " << mt1.molName << "\n";

    // Constructor with center-of-mass coordinate.
    std::cerr << "[INFO] Creating MolTemplate with (comCoord, interfaceList) constructor.\n";
    Coord com(5.0, 6.0, 7.0);
    MolTemplate mt2(com, ifaceList);
    EXPECT_DOUBLE_EQ(mt2.comCoord.x, 5.0) << "comCoord.x should be 5.0";
    EXPECT_DOUBLE_EQ(mt2.comCoord.y, 6.0) << "comCoord.y should be 6.0";
    EXPECT_DOUBLE_EQ(mt2.comCoord.z, 7.0) << "comCoord.z should be 7.0";
    std::cerr << "[INFO] mt2.comCoord = (" << mt2.comCoord.x << ", "
              << mt2.comCoord.y << ", " << mt2.comCoord.z << ")\n";
}

// ---------------------------------------------------------------------------
// Test MolTemplate::find_relIndex_from_absIndex().
// This maps an absolute State index back to the relative interface index.
// Note: on failure the real implementation calls exit(1), so we only test
// valid look-ups here to keep the test process alive.
// ---------------------------------------------------------------------------
void test_find_relIndex_from_absIndex()
{
    print_test_banner("MolTemplate::find_relIndex_from_absIndex", "class_moltemplate_functions.cpp");

    Coord c(0.0, 0.0, 0.0);

    // Build interfaces, each carrying a single state with a known absolute index.
    std::vector<Interface> ifaceList;

    // Interface 0, relative index 0, state absolute index 100.
    {
        std::vector<Interface::State> st;
        st.emplace_back("i0~s", 's', 100);
        Interface iface("i0", st, c);
        iface.index = 0; // relative index
        ifaceList.push_back(iface);
    }
    // Interface 1, relative index 1, state absolute index 200.
    {
        std::vector<Interface::State> st;
        st.emplace_back("i1~s", 's', 200);
        Interface iface("i1", st, c);
        iface.index = 1; // relative index
        ifaceList.push_back(iface);
    }

    MolTemplate mt("multiIface", ifaceList);

    // Absolute index 100 belongs to the first interface -> relative index 0.
    std::cerr << "[INFO] Looking up relative index for absolute state index 100.\n";
    int rel0 = mt.find_relIndex_from_absIndex(100);
    EXPECT_EQ(rel0, 0) << "Absolute index 100 should map to relative index 0";
    std::cerr << "[INFO] Result relative index = " << rel0 << "\n";

    // Absolute index 200 belongs to the second interface -> relative index 1.
    std::cerr << "[INFO] Looking up relative index for absolute state index 200.\n";
    int rel1 = mt.find_relIndex_from_absIndex(200);
    EXPECT_EQ(rel1, 1) << "Absolute index 200 should map to relative index 1";
    std::cerr << "[INFO] Result relative index = " << rel1 << "\n";
}

// ---------------------------------------------------------------------------
// Test MolTemplate::find_absIndex_from_relIndex().
// This maps (relative interface index, state identity) to the absolute
// State index. Again we only test valid look-ups to avoid exit(1).
// ---------------------------------------------------------------------------
void test_find_absIndex_from_relIndex()
{
    print_test_banner("MolTemplate::find_absIndex_from_relIndex", "class_moltemplate_functions.cpp");

    Coord c(0.0, 0.0, 0.0);
    std::vector<Interface> ifaceList;

    // Interface with relative index 0 and two states 'A' (abs 10) and 'B' (abs 11).
    {
        std::vector<Interface::State> st;
        st.emplace_back("phos~A", 'A', 10);
        st.emplace_back("phos~B", 'B', 11);
        Interface iface("phos", st, c);
        iface.index = 0;
        ifaceList.push_back(iface);
    }

    MolTemplate mt("stateMol", ifaceList);

    // State 'A' on relative interface 0 has absolute index 10.
    std::cerr << "[INFO] Looking up absolute index for state 'A' on interface 0.\n";
    int absA = mt.find_absIndex_from_relIndex(0, 'A');
    EXPECT_EQ(absA, 10) << "State 'A' on interface 0 should have absolute index 10";
    std::cerr << "[INFO] Result absolute index = " << absA << "\n";

    // State 'B' on relative interface 0 has absolute index 11.
    std::cerr << "[INFO] Looking up absolute index for state 'B' on interface 0.\n";
    int absB = mt.find_absIndex_from_relIndex(0, 'B');
    EXPECT_EQ(absB, 11) << "State 'B' on interface 0 should have absolute index 11";
    std::cerr << "[INFO] Result absolute index = " << absB << "\n";
}

// ---------------------------------------------------------------------------
// Test MolTemplate::set_value().
// This parser helper interprets a text 'line' according to a MolKeyword and
// mutates the MolTemplate accordingly. We exercise several keyword cases.
// ---------------------------------------------------------------------------
void test_set_value()
{
    print_test_banner("MolTemplate::set_value", "class_moltemplate_functions.cpp");

    Coord c(0.0, 0.0, 0.0);
    std::vector<Interface> ifaceList;
    ifaceList.emplace_back("iface0", c);
    MolTemplate mt("initial", ifaceList);

    // Keyword 0 (name): set the molecule name.
    std::cerr << "[INFO] set_value with keyword index 0 (name) -> 'newName'.\n";
    std::string nameLine = "newName";
    mt.set_value(nameLine, static_cast<MolKeyword>(0));
    EXPECT_EQ(mt.molName, "newName") << "molName should be updated to 'newName'";

    // Keyword 1 (copies): set the copy number.
    std::cerr << "[INFO] set_value with keyword index 1 (copies) -> 42.\n";
    std::string copiesLine = "42";
    mt.set_value(copiesLine, static_cast<MolKeyword>(1));
    EXPECT_EQ(mt.copies, 42) << "copies should be 42";

    // Keyword 2 (isRod): set boolean rod flag.
    std::cerr << "[INFO] set_value with keyword index 2 (isRod) -> true.\n";
    std::string rodLine = "true";
    mt.set_value(rodLine, static_cast<MolKeyword>(2));
    EXPECT_TRUE(mt.isRod) << "isRod should be true";

    // Keyword 8 (mass): set the mass value.
    std::cerr << "[INFO] set_value with keyword index 8 (mass) -> 3.5.\n";
    std::string massLine = "3.5";
    mt.set_value(massLine, static_cast<MolKeyword>(8));
    EXPECT_DOUBLE_EQ(mt.mass, 3.5) << "mass should be 3.5";

    // Keyword 12 (isImplicitLipid): when true, isLipid should also be set true.
    std::cerr << "[INFO] set_value with keyword index 12 (isImplicitLipid) -> true.\n";
    std::string implicitLine = "true";
    mt.set_value(implicitLine, static_cast<MolKeyword>(12));
    EXPECT_TRUE(mt.isImplicitLipid) << "isImplicitLipid should be true";
    EXPECT_TRUE(mt.isLipid) << "isLipid should follow isImplicitLipid";

    // Keyword 15 (outsideCompartment): when true, crossesCompartment true.
    std::cerr << "[INFO] set_value with keyword index 15 (outsideCompartment) -> true.\n";
    std::string outsideLine = "true";
    mt.set_value(outsideLine, static_cast<MolKeyword>(15));
    EXPECT_TRUE(mt.outsideCompartment) << "outsideCompartment should be true";
    EXPECT_TRUE(mt.crossesCompartment) << "crossesCompartment should follow outsideCompartment";

    // Keyword 17 (isPromoter): set boolean promoter flag.
    std::cerr << "[INFO] set_value with keyword index 17 (isPromoter) -> true.\n";
    std::string promoterLine = "true";
    mt.set_value(promoterLine, static_cast<MolKeyword>(17));
    EXPECT_TRUE(mt.isPromoter) << "isPromoter should be true";
}

// ---------------------------------------------------------------------------
// Test the free function find_molTypeIndex_from_ifaceIndex().
// It searches a list of MolTemplates for a state matching the given absolute
// interface index, and returns that template's molTypeIndex.
// Only valid look-ups are tested (invalid triggers exit(1)).
// ---------------------------------------------------------------------------
void test_find_molTypeIndex_from_ifaceIndex()
{
    print_test_banner("find_molTypeIndex_from_ifaceIndex", "class_moltemplate_functions.cpp");

    Coord c(0.0, 0.0, 0.0);

    // Build two MolTemplates each with a single interface/state.
    std::vector<MolTemplate> molTemplateList;

    // Template A: molTypeIndex 0, state absolute index 500.
    {
        std::vector<Interface::State> st;
        st.emplace_back("a~s", 's', 500);
        Interface iface("a", st, c);
        std::vector<Interface> ifaceList { iface };
        MolTemplate mtA("A", ifaceList);
        mtA.molTypeIndex = 0;
        molTemplateList.push_back(mtA);
    }
    // Template B: molTypeIndex 1, state absolute index 600.
    {
        std::vector<Interface::State> st;
        st.emplace_back("b~s", 's', 600);
        Interface iface("b", st, c);
        std::vector<Interface> ifaceList { iface };
        MolTemplate mtB("B", ifaceList);
        mtB.molTypeIndex = 1;
        molTemplateList.push_back(mtB);
    }

    // Absolute interface index 500 belongs to template A (molTypeIndex 0).
    std::cerr << "[INFO] Looking up molTypeIndex for absolute iface index 500.\n";
    int target500 = 500;
    size_t typeA = find_molTypeIndex_from_ifaceIndex(target500, molTemplateList);
    EXPECT_EQ(typeA, 0u) << "Iface index 500 should belong to molTypeIndex 0";
    std::cerr << "[INFO] Result molTypeIndex = " << typeA << "\n";

    // Absolute interface index 600 belongs to template B (molTypeIndex 1).
    std::cerr << "[INFO] Looking up molTypeIndex for absolute iface index 600.\n";
    int target600 = 600;
    size_t typeB = find_molTypeIndex_from_ifaceIndex(target600, molTemplateList);
    EXPECT_EQ(typeB, 1u) << "Iface index 600 should belong to molTypeIndex 1";
    std::cerr << "[INFO] Result molTypeIndex = " << typeB << "\n";
}

// ---------------------------------------------------------------------------
// Test MolTemplate::display() overloads.
// These functions only print to std::cout; there is no return value to check.
// We simply invoke them to ensure they execute without crashing and to
// demonstrate the produced output.
// ---------------------------------------------------------------------------
void mol_template_test_display_functions()
{
    print_test_banner("MolTemplate::display / display(name)", "class_moltemplate_functions.cpp");

    Coord c(1.0, 1.0, 1.0);

    // Build an interface with a single state.
    std::vector<Interface::State> singleState;
    singleState.emplace_back("i~s", 's', 0);
    Interface singleIface("iSingle", singleState, c);
    singleIface.index = 0;

    // Build an interface with multiple states.
    std::vector<Interface::State> multiStates;
    multiStates.emplace_back("m~A", 'A', 1);
    multiStates.emplace_back("m~B", 'B', 2);
    Interface multiIface("iMulti", multiStates, c);
    multiIface.index = 1;

    std::vector<Interface> ifaceList { singleIface, multiIface };
    MolTemplate mt("displayMol", ifaceList);

    // Invoke display() - if it returns we consider it a success.
    std::cerr << "[INFO] Calling MolTemplate::display() (output follows below).\n";
    mt.display();
    SUCCEED() << "MolTemplate::display() executed without crashing";

    // Invoke display(name) overload.
    std::cerr << "[INFO] Calling MolTemplate::display(name) (output follows below).\n";
    mt.display("displayMolByName");
    SUCCEED() << "MolTemplate::display(name) executed without crashing";
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers. Each TEST simply calls one of the test_* functions.
// Using non-fatal EXPECT_* assertions inside those functions ensures all
// checks run even if some fail.
// ---------------------------------------------------------------------------

TEST(MolTemplateTest, StateConstructors) { test_State_constructors(); }
TEST(MolTemplateTest, InterfaceConstructorsAndNaming) { test_Interface_constructors_and_naming(); }
TEST(MolTemplateTest, MolTemplateConstructors) { test_MolTemplate_constructors(); }
TEST(MolTemplateTest, FindRelIndexFromAbsIndex) { test_find_relIndex_from_absIndex(); }
TEST(MolTemplateTest, FindAbsIndexFromRelIndex) { test_find_absIndex_from_relIndex(); }
TEST(MolTemplateTest, SetValue) { test_set_value(); }
TEST(MolTemplateTest, FindMolTypeIndexFromIfaceIndex) { test_find_molTypeIndex_from_ifaceIndex(); }
TEST(MolTemplateTest, DisplayFunctions) { mol_template_test_display_functions(); }