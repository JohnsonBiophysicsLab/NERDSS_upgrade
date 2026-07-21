
// =====================================================================================
// Unit tests for cluster / trajectory functions defined in the provided .cpp source.
//
// The source file (class_Cluster.cpp style) defines:
//   - ClusterPair::ClusterPair()             (default constructor)
//   - ClusterPair::ClusterPair(int, int)     (parameterized constructor)
//   - cluster_one_complex(...)               (find cross-partner pairs for one complex)
//   - define_cluster_pairs(...)              (descend through partners collecting pairs)
//   - resample_traj(...)                     (resample trajectory of complexes in pairs)
//
// These tests focus on the behavior we can exercise in isolation:
//   * The ClusterPair constructors.
//   * cluster_one_complex / define_cluster_pairs pair-building logic.
//
// We build minimal Molecule/Complex/ForwardRxn/Membrane objects, invoke the
// functions, and verify the resulting pairList contents.
//
// The test framework is GoogleTest, but we additionally use small helper
// functions (require_close / require_true) that print verbose diagnostics to
// stderr and call std::exit(1) on failure, as requested.
// =====================================================================================

#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

// Headers describing the types used by the source under test.
#include "classes/class_Cluster.hpp"
#include "boundary_conditions/reflect_functions.hpp"
#include "trajectory_functions/trajectory_functions.hpp"

// -------------------------------------------------------------------------------------
// Forward declarations of the functions under test (as defined in the source file).
// -------------------------------------------------------------------------------------
void cluster_one_complex(int k1, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    std::vector<ClusterPair>& pairList, std::vector<int>& finished);

void define_cluster_pairs(int p1, std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList, const std::vector<ForwardRxn>& forwardRxns,
    std::vector<ClusterPair>& pairList);

// =====================================================================================
// Verbose helper assertion functions.
// These print what they check and, on failure, print an error and abort the program.
// =====================================================================================

// Compare two doubles within a small tolerance; abort on failure.
static void require_close(double actual, double expected, const std::string& label)
{
    const double tol = 1e-9;
    std::cerr << "    [require_close] " << label
              << " : actual=" << actual << " expected=" << expected << std::endl;
    if (std::fabs(actual - expected) > tol) {
        std::cerr << "    FAILED: " << label << " differ by "
                  << std::fabs(actual - expected) << " (tol=" << tol << ")" << std::endl;
        std::exit(1);
    }
}

// Verify a boolean condition; abort on failure.
static void require_true(bool condition, const std::string& label)
{
    std::cerr << "    [require_true]  " << label
              << " : " << (condition ? "true" : "false") << std::endl;
    if (!condition) {
        std::cerr << "    FAILED: condition '" << label << "' was not true" << std::endl;
        std::exit(1);
    }
}

// =====================================================================================
// Helpers to build a minimal simulation setup for the pair-building functions.
// =====================================================================================

// Build a very simple two-molecule / two-complex scenario where molecule 0
// (in complex 0) cross-reacts with molecule 1 (in complex 1) via reaction 0.
// This is the minimal input needed to exercise cluster_one_complex.
static void build_two_molecule_scenario(
    std::vector<Molecule>& moleculeList,
    std::vector<Complex>& complexList,
    std::vector<ForwardRxn>& forwardRxns)
{
    // --- Set up the single forward reaction (index 0). ---
    // Interface index 0 of reactant 0 binds to interface index 0 of reactant 1.
    ForwardRxn rxn {};
    rxn.bindRadius = 1.5;
    // reactantListNew must contain 2 entries so the code can pick the partner interface.
    rxn.reactantListNew.resize(2);
    rxn.reactantListNew[0].relIfaceIndex = 0; // interface on molecule p1
    rxn.reactantListNew[1].relIfaceIndex = 0; // interface on molecule p2
    forwardRxns.clear();
    forwardRxns.push_back(rxn);

    // --- Set up two molecules. ---
    moleculeList.clear();
    moleculeList.resize(2);

    // Molecule 0 lives in complex 0 and crosses with molecule 1 through reaction 0.
    moleculeList[0].myComIndex = 0;
    moleculeList[0].trajStatus = TrajStatus::none;
    moleculeList[0].crossbase = { 1 };                 // partner molecule index
    moleculeList[0].mycrossint = { 0 };                // reacting interface index on mol 0
    moleculeList[0].crossrxn = { std::array<int, 3> { 0, 0, 0 } }; // rxn index in [0]

    // Molecule 1 lives in complex 1; it has no outgoing cross list of its own here.
    moleculeList[1].myComIndex = 1;
    moleculeList[1].trajStatus = TrajStatus::none;
    moleculeList[1].crossbase = {};                    // no cross partners recorded
    moleculeList[1].mycrossint = {};
    moleculeList[1].crossrxn = {};

    // --- Set up the two complexes, each containing one molecule. ---
    complexList.clear();
    complexList.resize(2);

    complexList[0].memberList = { 0 };
    complexList[0].OnSurface = false; // volume complex -> high priority expected

    complexList[1].memberList = { 1 };
    complexList[1].OnSurface = false;
}

// =====================================================================================
// Test group 1: ClusterPair constructors.
// =====================================================================================
void test_clusterpair_constructors()
{
    std::cerr << "\n[TEST] test_clusterpair_constructors "
              << "(file: class_Cluster.cpp, ClusterPair constructors)" << std::endl;

    // Default constructor should set p1 = -1 (sentinel meaning "unset").
    std::cerr << "  Checking default constructor sets p1 = -1..." << std::endl;
    ClusterPair defaultPair;
    require_close(static_cast<double>(defaultPair.p1), -1.0,
        "default ClusterPair p1 sentinel");

    // Parameterized constructor should store both provided indices.
    std::cerr << "  Checking parameterized constructor stores p1 and p2..." << std::endl;
    ClusterPair setPair(7, 42);
    require_close(static_cast<double>(setPair.p1), 7.0, "parameterized ClusterPair p1");
    require_close(static_cast<double>(setPair.p2), 42.0, "parameterized ClusterPair p2");

    std::cerr << "  PASSED: ClusterPair constructors behave as expected." << std::endl;
}

// =====================================================================================
// Test group 2: cluster_one_complex builds the correct cross pairs.
// =====================================================================================
void test_cluster_one_complex()
{
    std::cerr << "\n[TEST] test_cluster_one_complex "
              << "(file: class_Cluster.cpp, cluster_one_complex)" << std::endl;

    // Build the minimal two-molecule / two-complex scenario.
    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
    build_two_molecule_scenario(moleculeList, complexList, forwardRxns);

    std::vector<ClusterPair> pairList;
    std::vector<int> finished;

    // Run for complex 0 (which contains molecule 0 that crosses with molecule 1).
    std::cerr << "  Calling cluster_one_complex(k1=0, ...)..." << std::endl;
    cluster_one_complex(0, moleculeList, complexList, forwardRxns, pairList, finished);

    // We expect exactly one cross-complex pair to be recorded.
    std::cerr << "  pairList.size() = " << pairList.size() << std::endl;
    require_true(pairList.size() == 1, "cluster_one_complex produced exactly one pair");

    // Verify the details of the discovered pair.
    const ClusterPair& pr = pairList[0];
    require_true(pr.p1 == 0, "pair.p1 == 0 (source molecule)");
    require_true(pr.p2 == 1, "pair.p2 == 1 (partner molecule)");
    require_true(pr.k1 == 0, "pair.k1 == 0 (source complex)");
    require_true(pr.k2 == 1, "pair.k2 == 1 (partner complex)");
    require_true(pr.i1 == 0, "pair.i1 == 0 (reacting interface on p1)");
    require_true(pr.i2 == 0, "pair.i2 == 0 (partner interface on p2)");
    require_close(pr.bindrad, forwardRxns[0].bindRadius, "pair.bindrad == rxn bindRadius");

    // Both complexes are volume complexes -> priority should be high (1), memtest 0.
    require_true(pr.priority == 1, "pair.priority == 1 (both volume complexes)");
    require_true(pr.memtest == 0, "pair.memtest == 0 (neither complex on surface)");

    // The processed complex (k1=0) must be recorded in finished.
    require_true(!finished.empty() && finished[0] == 0,
        "finished contains the processed complex index 0");

    std::cerr << "  PASSED: cluster_one_complex built the expected cross pair." << std::endl;
}

// =====================================================================================
// Test group 3: cluster_one_complex marks surface complexes with reduced priority.
// =====================================================================================
void test_cluster_one_complex_surface_priority()
{
    std::cerr << "\n[TEST] test_cluster_one_complex_surface_priority "
              << "(file: class_Cluster.cpp, cluster_one_complex surface logic)" << std::endl;

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
    build_two_molecule_scenario(moleculeList, complexList, forwardRxns);

    // Mark BOTH complexes as being on the surface.
    std::cerr << "  Marking both complexes as OnSurface = true..." << std::endl;
    complexList[0].OnSurface = true;
    complexList[1].OnSurface = true;

    std::vector<ClusterPair> pairList;
    std::vector<int> finished;

    cluster_one_complex(0, moleculeList, complexList, forwardRxns, pairList, finished);

    require_true(pairList.size() == 1, "surface scenario still produces one pair");
    const ClusterPair& pr = pairList[0];

    // On-surface complexes -> low priority (0) and memtest == 1 when both are on surface.
    require_true(pr.priority == 0, "pair.priority == 0 (a complex is on surface)");
    require_true(pr.memtest == 1, "pair.memtest == 1 (both complexes on surface)");

    std::cerr << "  PASSED: surface priority / memtest flags set correctly." << std::endl;
}

// =====================================================================================
// Test group 4: define_cluster_pairs descends through complexes without duplicating.
// =====================================================================================
void test_define_cluster_pairs()
{
    std::cerr << "\n[TEST] test_define_cluster_pairs "
              << "(file: class_Cluster.cpp, define_cluster_pairs)" << std::endl;

    std::vector<Molecule> moleculeList;
    std::vector<Complex> complexList;
    std::vector<ForwardRxn> forwardRxns;
    build_two_molecule_scenario(moleculeList, complexList, forwardRxns);

    std::vector<ClusterPair> pairList;

    // Starting from molecule 0, define_cluster_pairs should collect the same single
    // cross pair (0->1) and then follow into complex 1, which has no cross partners,
    // so no additional pairs should be added.
    std::cerr << "  Calling define_cluster_pairs(p1=0, ...)..." << std::endl;
    define_cluster_pairs(0, moleculeList, complexList, forwardRxns, pairList);

    std::cerr << "  pairList.size() = " << pairList.size() << std::endl;
    require_true(pairList.size() == 1,
        "define_cluster_pairs collects exactly one pair (no duplicates)");

    const ClusterPair& pr = pairList[0];
    require_true(pr.p1 == 0 && pr.p2 == 1,
        "define_cluster_pairs pair connects molecule 0 and molecule 1");

    std::cerr << "  PASSED: define_cluster_pairs traversed complexes correctly." << std::endl;
}

// =====================================================================================
// GoogleTest wrappers. Each simply invokes the verbose test_* function above.
// (If any require_* helper fails it calls std::exit(1), aborting the whole run.)
// =====================================================================================

TEST(ClusterTests, ClusterPairConstructors)          { test_clusterpair_constructors(); }
TEST(ClusterTests, ClusterOneComplex)                { test_cluster_one_complex(); }
TEST(ClusterTests, ClusterOneComplexSurfacePriority) { test_cluster_one_complex_surface_priority(); }
TEST(ClusterTests, DefineClusterPairs)               { test_define_cluster_pairs(); }

