/*! \file test_class_rxn_functions.cpp
 *
 * ### Purpose
 * ***
 * Unit tests for the reaction class functions defined in class_rxn_functions.cpp
 * (part of class_Rxns.hpp). These tests exercise:
 *   - The RxnIface constructor and its operator<< / operator== overloads
 *   - The ForwardRxn::Angles::display() method
 *   - The static counters maintained by RxnBase (numberOfRxns / totRxnSpecies)
 *
 * NOTE:
 *   The full constructors for ForwardRxn / BackRxn / CreateDestructRxn require a
 *   fully-parsed ParsedRxn along with a valid MolTemplate list, which is a large
 *   dependency. Therefore these tests focus on the smaller, self-contained pieces
 *   (RxnIface behavior and static bookkeeping) that can be tested in isolation.
 */

#include "classes/class_Rxns.hpp"          // RxnBase, RxnIface, ForwardRxn, etc.
#include "classes/class_bngl_parser.hpp"   // ParsedRxn (needed transitively)

#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>

// ---------------------------------------------------------------------------
// Helper functions for verbose reporting of test criteria and results.
// ---------------------------------------------------------------------------

/*! \brief Check that two floating point values are close, printing progress.
 *  Prints the label, expected, and actual values to stderr.*/
static void require_close(double actual, double expected, const std::string& label)
{
    // Tolerance used for floating-point comparison.
    const double tol = 1e-9;
    std::cerr << "  [CHECK] " << label
              << " : expected = " << expected
              << ", actual = " << actual << std::endl;
    if (std::fabs(actual - expected) > tol) {
        std::cerr << "  [FAIL ] " << label
                  << " : difference " << std::fabs(actual - expected)
                  << " exceeds tolerance " << tol << std::endl;
    
    }
    EXPECT_FALSE(std::fabs(actual - expected) > tol);
    std::cerr << "  [PASS ] " << label << std::endl;
}

/*! \brief Check that a boolean condition holds, printing progress.
 *  Prints the label to stderr. */
static void require_true(bool condition, const std::string& label)
{
    std::cerr << "  [CHECK] " << label << std::endl;
    if (!condition) {
        std::cerr << "  [FAIL ] " << label << " : condition was false" << std::endl;
        
    }
    EXPECT_FALSE(!condition);
    std::cerr << "  [PASS ] " << label << std::endl;
}

// ---------------------------------------------------------------------------
// Tests for the RxnIface class (constructor + member storage)
// ---------------------------------------------------------------------------

/*! \brief Verify the RxnIface constructor stores all the supplied fields. */
void test_rxniface_constructor()
{
    std::cerr << "\n=== Testing RxnIface constructor "
                 "(class_rxn_functions.cpp) ===" << std::endl;

    // Build an interface: named "site", molTypeIndex 2, absIndex 5, relIndex 1,
    // requires state 'P' (phosphorylated), and requires a bound interaction.
    RxnIface iface("site", 2, 5, 1, 'P', true);

    // Confirm every field is stored exactly as passed to the constructor.
    require_true(iface.ifaceName == "site", "ifaceName stored correctly");
    require_true(iface.molTypeIndex == 2, "molTypeIndex stored correctly");
    require_true(iface.absIfaceIndex == 5, "absIfaceIndex stored correctly");
    require_true(iface.relIfaceIndex == 1, "relIfaceIndex stored correctly");
    require_true(iface.requiresState == 'P', "requiresState stored correctly");
    require_true(iface.requiresInteraction == true,
                 "requiresInteraction stored correctly");
}

// ---------------------------------------------------------------------------
// Tests for the RxnIface operator<< overload
// ---------------------------------------------------------------------------

/*! \brief Verify the streaming operator produces the expected human-readable text. */
void test_rxniface_ostream()
{
    std::cerr << "\n=== Testing RxnIface operator<< "
                 "(class_rxn_functions.cpp) ===" << std::endl;

    // Case 1: interface that requires a state and a bound interaction.
    {
        RxnIface iface("A", 0, 3, 0, 'U', true);
        std::ostringstream oss;
        oss << iface;  // Invoke the custom operator<<.
        std::string out = oss.str();
        std::cerr << "  produced string: " << out << std::endl;

        // The output must contain the interface name and absolute index.
        require_true(out.find("iface name: A") != std::string::npos,
                     "output contains interface name");
        require_true(out.find("iface index: 3") != std::string::npos,
                     "output contains absolute index");
        // Because requiresState != '\0', the state should be printed.
        require_true(out.find("requires state U") != std::string::npos,
                     "output contains required state");
        // Because requiresInteraction is true, "requires interaction" is printed.
        require_true(out.find("requires interaction") != std::string::npos,
                     "output notes required interaction");
    }

    // Case 2: interface that is free (no interaction) and no state requirement.
    {
        // '\0' state means no state requirement is printed; false means "required Free".
        RxnIface iface("B", 1, 7, 2, '\0', false);
        std::ostringstream oss;
        oss << iface;
        std::string out = oss.str();
        std::cerr << "  produced string: " << out << std::endl;

        require_true(out.find("iface name: B") != std::string::npos,
                     "output contains interface name (free case)");
        // No state requirement should appear.
        require_true(out.find("requires state") == std::string::npos,
                     "output omits state when none required");
        // Because requiresInteraction is false, "required Free" is printed.
        require_true(out.find("required Free") != std::string::npos,
                     "output notes free interface");
    }
}

// ---------------------------------------------------------------------------
// Tests for the RxnIface::operator== (comparison against a Molecule::Iface)
// ---------------------------------------------------------------------------

/*! \brief Verify equality logic between a RxnIface and a Molecule::Iface. */
void test_rxniface_equality()
{
    std::cerr << "\n=== Testing RxnIface::operator==(Molecule::Iface) "
                 "(class_rxn_functions.cpp) ===" << std::endl;

    // Construct a reaction interface that we will compare against molecule ifaces.
    // molTypeIndex = 1, relIndex = 0, requiresState = 'S', requiresInteraction = false.
    RxnIface rxnIface("site", 1, 4, 0, 'S', false);

    // Build a molecule interface that matches on all compared fields.
    Molecule::Iface matching{};
    matching.molTypeIndex = 1;    // matches rxnIface.molTypeIndex
    matching.relIndex     = 0;    // matches rxnIface.relIfaceIndex
    matching.isBound      = false; // matches rxnIface.requiresInteraction
    matching.stateIden    = 'S';  // matches rxnIface.requiresState

    // The comparison should succeed for a fully matching interface.
    require_true(rxnIface == matching, "matching interface compares equal");

    // Now build several mismatching interfaces; each should compare NOT equal.

    // Mismatch on molTypeIndex.
    {
        Molecule::Iface m = matching;
        m.molTypeIndex = 99;
        require_true(!(rxnIface == m), "different molTypeIndex compares unequal");
    }

    // Mismatch on relIndex.
    {
        Molecule::Iface m = matching;
        m.relIndex = 5;
        require_true(!(rxnIface == m), "different relIndex compares unequal");
    }

    // Mismatch on isBound (interaction requirement).
    {
        Molecule::Iface m = matching;
        m.isBound = true;
        require_true(!(rxnIface == m), "different isBound compares unequal");
    }

    // Mismatch on state identifier.
    {
        Molecule::Iface m = matching;
        m.stateIden = 'X';
        require_true(!(rxnIface == m), "different stateIden compares unequal");
    }
}

// ---------------------------------------------------------------------------
// Tests for the ForwardRxn::Angles::display() method
// ---------------------------------------------------------------------------

/*! \brief Verify that setting angle members and calling display() works and
 *         stores the values as expected. */
void test_forwardrxn_angles()
{
    std::cerr << "\n=== Testing ForwardRxn::Angles::display() "
                 "(class_rxn_functions.cpp) ===" << std::endl;

    // Populate an Angles struct with distinct known values.
    ForwardRxn::Angles angles{};
    angles.theta1 = 1.1;
    angles.theta2 = 2.2;
    angles.phi1   = 3.3;
    angles.phi2   = 4.4;
    angles.omega  = 5.5;

    // Confirm the stored values are exactly what we set.
    require_close(angles.theta1, 1.1, "theta1 stored");
    require_close(angles.theta2, 2.2, "theta2 stored");
    require_close(angles.phi1, 3.3, "phi1 stored");
    require_close(angles.phi2, 4.4, "phi2 stored");
    require_close(angles.omega, 5.5, "omega stored");

    // Exercise the display() method to ensure it runs without crashing.
    std::cerr << "  invoking Angles::display() (output goes to stdout):" << std::endl;
    angles.display();
    require_true(true, "Angles::display() executed without error");
}

// ---------------------------------------------------------------------------
// Tests for the static RxnBase bookkeeping counters.
// ---------------------------------------------------------------------------

/*! \brief Verify that RxnBase static members are accessible and mutable, and
 *         behave like a shared counter across the reaction hierarchy. */
void test_rxnbase_static_counters()
{
    std::cerr << "\n=== Testing RxnBase static counters "
                 "(class_rxn_functions.cpp) ===" << std::endl;

    // Snapshot the current value of the shared reaction counter.
    unsigned before = RxnBase::numberOfRxns;
    std::cerr << "  numberOfRxns before increment: " << before << std::endl;

    // Manually increment and verify the change is observed.
    ++RxnBase::numberOfRxns;
    require_true(RxnBase::numberOfRxns == before + 1,
                 "numberOfRxns increments correctly");

    // Restore the counter so we do not perturb other tests in the suite.
    RxnBase::numberOfRxns = before;
    require_true(RxnBase::numberOfRxns == before,
                 "numberOfRxns restored to original value");

    // Verify the species counter is likewise accessible.
    int speciesBefore = RxnBase::totRxnSpecies;
    RxnBase::totRxnSpecies = speciesBefore + 3;
    require_true(RxnBase::totRxnSpecies == speciesBefore + 3,
                 "totRxnSpecies is mutable");
    RxnBase::totRxnSpecies = speciesBefore; // restore
    require_true(RxnBase::totRxnSpecies == speciesBefore,
                 "totRxnSpecies restored to original value");
}

// ---------------------------------------------------------------------------
// GoogleTest wrappers: each TEST drives one of the verbose test_* functions.
// ---------------------------------------------------------------------------

// The RxnIface constructor stores all constructor arguments.
TEST(ClassRxnFunctions, RxnIfaceConstructor) { test_rxniface_constructor(); }

// The RxnIface operator<< produces a correct human-readable description.
TEST(ClassRxnFunctions, RxnIfaceOstream) { test_rxniface_ostream(); }

// The RxnIface operator==(Molecule::Iface) compares the correct fields.
TEST(ClassRxnFunctions, RxnIfaceEquality) { test_rxniface_equality(); }

// The ForwardRxn::Angles struct stores/displays its values.
TEST(ClassRxnFunctions, ForwardRxnAngles) { test_forwardrxn_angles(); }

// The static RxnBase counters behave as expected.
TEST(ClassRxnFunctions, RxnBaseStaticCounters) { test_rxnbase_static_counters(); }