
/*! \file test_class_Observable.cpp
 *
 * \brief Unit tests for the Observable class member functions defined in
 *        class_Observable.cpp.
 *
 * This test suite exercises:
 *   - Observable::Constituent::Iface::operator==(const RxnIface&)
 *   - Observable::operator==(const Molecule&)
 *
 * The tests print verbose information to stderr so it is clear which
 * function in which source file is being tested and what each test does.
 */

#include <gtest/gtest.h>
#include <iostream>

// Header that declares the Observable class and its nested types, plus
// the Molecule and RxnIface types referenced by the operators under test.
#include "classes/class_Observable.hpp"

using namespace SpeciesTracker;

/*
 * ------------------------------------------------------------------------
 * Tests for Observable::Constituent::Iface::operator==(const RxnIface&)
 * ------------------------------------------------------------------------
 *
 * The operator (from class_Observable.cpp) returns:
 *     (relIndex == rxnIface.relIfaceIndex) && (isBound == rxnIface.relIfaceIndex)
 *
 * Note: The second comparison compares 'isBound' (a bool) against
 * 'relIfaceIndex' (an integer). This is exactly the behaviour we test,
 * since the tests should reflect the actual implementation.
 */

//! Test that the Iface equality operator returns true when both conditions hold.
void test_Iface_equal_true()
{
    std::cerr << "[TEST] class_Observable.cpp :: "
                 "Observable::Constituent::Iface::operator==\n"
              << "       Verifying operator== returns TRUE when relIndex matches "
                 "relIfaceIndex and isBound matches relIfaceIndex.\n";

    // Build an Iface with relIndex == 1 and isBound == true (which converts to 1).
    Observable::Constituent::Iface iface;
    iface.relIndex = 1;
    iface.isBound  = true; // true == 1, so both comparisons should succeed.

    // Build a matching RxnIface with relIfaceIndex == 1.
    RxnIface rxnIface;
    rxnIface.relIfaceIndex = 1;

    // Both (relIndex == 1) and (isBound(true->1) == 1) are true.
    bool result = (iface == rxnIface);
    std::cerr << "       relIndex=1, isBound=true, relIfaceIndex=1 -> result="
              << std::boolalpha << result << "\n";

    // We expect equality to hold.
    EXPECT_TRUE(result)
        << "Expected Iface == RxnIface to be true when relIndex and "
           "isBound both equal relIfaceIndex.";
}

//! Test that the operator returns false when relIndex differs from relIfaceIndex.
void test_Iface_relIndex_mismatch()
{
    std::cerr << "[TEST] class_Observable.cpp :: "
                 "Observable::Constituent::Iface::operator==\n"
              << "       Verifying operator== returns FALSE when relIndex does "
                 "not match relIfaceIndex.\n";

    // relIndex (2) does not match relIfaceIndex (5), so operator== must be false.
    Observable::Constituent::Iface iface;
    iface.relIndex = 2;
    iface.isBound  = true;

    RxnIface rxnIface;
    rxnIface.relIfaceIndex = 5;

    bool result = (iface == rxnIface);
    std::cerr << "       relIndex=2, relIfaceIndex=5 -> result="
              << std::boolalpha << result << "\n";

    // Because relIndex != relIfaceIndex, the AND expression short-circuits false.
    EXPECT_FALSE(result)
        << "Expected Iface == RxnIface to be false when relIndex differs "
           "from relIfaceIndex.";
}

//! Test the second condition: isBound compared against relIfaceIndex.
void test_Iface_isBound_condition()
{
    std::cerr << "[TEST] class_Observable.cpp :: "
                 "Observable::Constituent::Iface::operator==\n"
              << "       Verifying the second comparison (isBound == "
                 "relIfaceIndex) affects the result.\n";

    // Here relIndex matches relIfaceIndex, but isBound (false -> 0) does NOT
    // equal relIfaceIndex (1), so the operator should return false.
    Observable::Constituent::Iface iface;
    iface.relIndex = 1;
    iface.isBound  = false; // false converts to 0, which != 1.

    RxnIface rxnIface;
    rxnIface.relIfaceIndex = 1;

    bool result = (iface == rxnIface);
    std::cerr << "       relIndex=1, isBound=false(0), relIfaceIndex=1 -> result="
              << std::boolalpha << result << "\n";

    // relIndex matches but isBound(0) != relIfaceIndex(1) so result is false.
    EXPECT_FALSE(result)
        << "Expected Iface == RxnIface to be false when isBound (0) does not "
           "equal relIfaceIndex (1).";
}

/*
 * ------------------------------------------------------------------------
 * Tests for Observable::operator==(const Molecule&)
 * ------------------------------------------------------------------------
 *
 * The operator (from class_Observable.cpp):
 *   1. Returns false if constituentList[0].molTypeIndex != mol.molTypeIndex.
 *   2. Otherwise iterates over constituentList[0].interfaceList and returns
 *      false if mol.interfaceList[obsIface.relIndex].index != obsIface.absIndex.
 *   3. Returns true if all checks pass.
 *
 * Helper functions build a matching or mismatching Observable/Molecule pair.
 */

//! Build a basic Observable with a single constituent and one interface.
static Observable makeObservable(int molTypeIndex, int relIndex, int absIndex)
{
    Observable obs;

    // Create a single constituent for the observable.
    Observable::Constituent constituent;
    constituent.molTypeIndex = molTypeIndex;

    // Add a single interface entry to the constituent.
    Observable::Constituent::Iface iface;
    iface.relIndex = relIndex;
    iface.absIndex = absIndex;
    constituent.interfaceList.push_back(iface);

    // Store the constituent so constituentList[0] is valid.
    obs.constituentList.push_back(constituent);

    return obs;
}

//! Build a Molecule with the given molTypeIndex and one interface entry.
static Molecule makeMolecule(int molTypeIndex, int ifaceIndexValue)
{
    Molecule mol;
    mol.molTypeIndex = molTypeIndex;

    // Add a single interface whose .index will be compared against absIndex.
    Molecule::Iface molIface;
    molIface.index = ifaceIndexValue;
    mol.interfaceList.push_back(molIface);

    return mol;
}

//! Test that Observable == Molecule is true for a fully matching pair.
void test_Observable_equal_true()
{
    std::cerr << "[TEST] class_Observable.cpp :: Observable::operator==\n"
              << "       Verifying operator== returns TRUE when molTypeIndex "
                 "and all interface indices match.\n";

    // Observable expects molecule type 3 and interface[0].index == 7.
    Observable obs = makeObservable(/*molTypeIndex=*/3, /*relIndex=*/0,
                                    /*absIndex=*/7);

    // Molecule of type 3 whose interface[0].index == 7 matches.
    Molecule mol = makeMolecule(/*molTypeIndex=*/3, /*ifaceIndexValue=*/7);

    bool result = (obs == mol);
    std::cerr << "       molTypeIndex match(3==3), absIndex(7)==mol.iface.index(7)"
                 " -> result=" << std::boolalpha << result << "\n";

    // All comparisons succeed, so we expect true.
    EXPECT_TRUE(result)
        << "Expected Observable == Molecule to be true when types and interface "
           "indices all match.";
}

//! Test that Observable == Molecule is false when molTypeIndex differs.
void test_Observable_molType_mismatch()
{
    std::cerr << "[TEST] class_Observable.cpp :: Observable::operator==\n"
              << "       Verifying operator== returns FALSE when molTypeIndex "
                 "does not match.\n";

    // Observable expects molecule type 3.
    Observable obs = makeObservable(/*molTypeIndex=*/3, /*relIndex=*/0,
                                    /*absIndex=*/7);

    // Molecule is type 9, which does not match the observable's type 3.
    Molecule mol = makeMolecule(/*molTypeIndex=*/9, /*ifaceIndexValue=*/7);

    bool result = (obs == mol);
    std::cerr << "       molTypeIndex mismatch (3 != 9) -> result="
              << std::boolalpha << result << "\n";

    // Early return false because molTypeIndex differs.
    EXPECT_FALSE(result)
        << "Expected Observable == Molecule to be false when molTypeIndex "
           "differs.";
}

//! Test that Observable == Molecule is false when interface index differs.
void test_Observable_iface_mismatch()
{
    std::cerr << "[TEST] class_Observable.cpp :: Observable::operator==\n"
              << "       Verifying operator== returns FALSE when an interface "
                 "index does not match.\n";

    // Observable expects molecule type 3 and interface[0].index == 7.
    Observable obs = makeObservable(/*molTypeIndex=*/3, /*relIndex=*/0,
                                    /*absIndex=*/7);

    // Molecule is type 3 (matches) but interface[0].index == 42 (does not match 7).
    Molecule mol = makeMolecule(/*molTypeIndex=*/3, /*ifaceIndexValue=*/42);

    bool result = (obs == mol);
    std::cerr << "       molType match(3), but absIndex(7) != mol.iface.index(42)"
                 " -> result=" << std::boolalpha << result << "\n";

    // Type matches but interface index differs, so result must be false.
    EXPECT_FALSE(result)
        << "Expected Observable == Molecule to be false when interface index "
           "differs.";
}

/*
 * ------------------------------------------------------------------------
 * GoogleTest wrappers - each TEST invokes one of the named test_* helpers.
 * Using EXPECT_* inside the helpers keeps failures non-fatal so that all
 * tests run even if some fail.
 * ------------------------------------------------------------------------
 */

// Iface::operator== tests.
TEST(ObservableIfaceTest, EqualReturnsTrue)      { test_Iface_equal_true(); }
TEST(ObservableIfaceTest, RelIndexMismatch)      { test_Iface_relIndex_mismatch(); }
TEST(ObservableIfaceTest, IsBoundCondition)      { test_Iface_isBound_condition(); }

// Observable::operator== tests.
TEST(ObservableTest, EqualReturnsTrue)           { test_Observable_equal_true(); }
TEST(ObservableTest, MolTypeMismatch)            { test_Observable_molType_mismatch(); }
TEST(ObservableTest, IfaceMismatch)              { test_Observable_iface_mismatch(); }