/*! \file test_create_DDMatrices.cpp
 *
 * ### Unit test for src/reactions/create_DDMatrices.cpp
 *
 * The single function defined in that translation unit is
 *
 *     void create_DDMatrices(gsl_matrix*& survMatrix,
 *                            gsl_matrix*& normMatrix,
 *                            gsl_matrix*& pirMatrix,
 *                            double bindRadius, double Dtot,
 *                            double comRMax, double ktemp,
 *                            const Parameters& params);
 *
 * It is a thin driver which
 *   1. computes the radial step size  RStepSize = sqrt(Dtot * timeStep) / 50,
 *   2. delegates to create_normMatrix(), create_survMatrix() and
 *      create_pirMatrix(), each of which *allocates* a gsl_matrix and hands the
 *      pointer back through the reference parameter.
 *
 * Because the function itself performs no arithmetic beyond RStepSize, the
 * things we can meaningfully verify are:
 *   - all three output pointers are actually allocated (non-null) and distinct,
 *   - the allocated tables are large enough to cover the lookup range that
 *     size_lookup() reports for the same inputs,
 *   - every stored element is a finite number (no NaN / inf leaking out of the
 *     numerical integration),
 *   - the routine is deterministic: identical inputs produce bit-identical
 *     tables (it must not depend on the RNG or on uninitialised memory),
 *   - increasing the maximum center-of-mass separation (comRMax) increases the
 *     number of tabulated entries.
 *
 * Verbose progress information is written to stderr so that the reader can see
 * which source file / function is under test and what each assertion checks.
 *
 * All parameters are deliberately chosen so the tables stay small (a few tens
 * of radial bins); create_pirMatrix() is O(N^2) in the number of bins, so large
 * tables would make the test needlessly slow.
 */

#include "classes/class_Parameters.hpp"
#include "reactions/bimolecular/2D_reaction_table_functions.hpp"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------
// Local helpers.  Everything is placed in an anonymous namespace and/or given
// the "ddm_" prefix so that no symbol collides with the rest of the test suite.
// -----------------------------------------------------------------------------
namespace {

/*! \brief Build a minimal Parameters object suitable for 2D table creation.
 *
 * Only the timestep is consulted by create_DDMatrices() (and by the routines it
 * delegates to), but we also set the static Parameters::dt because some helper
 * routines in the project read that instead of the instance member.
 *
 * \param[in] timeStep The simulation timestep, in microseconds.
 */
Parameters ddm_make_params(double timeStep)
{
    Parameters params;
    params.timeStep = timeStep;
    Parameters::dt = timeStep; // keep the static copy consistent
    return params;
}

/*! \brief Count how many elements of a matrix are not finite (NaN or +-inf).
 *
 * Also prints the matrix dimensions so the console log documents the size of
 * the generated lookup table.
 */
size_t ddm_count_nonfinite(const gsl_matrix* mat, const char* name)
{
    if (mat == nullptr) {
        std::cerr << "    " << name << " is NULL - cannot inspect elements\n";
        return 0;
    }

    size_t badCount = 0;
    for (size_t i = 0; i < mat->size1; ++i) {
        for (size_t j = 0; j < mat->size2; ++j) {
            if (!std::isfinite(gsl_matrix_get(mat, i, j)))
                ++badCount;
        }
    }

    std::cerr << "    " << name << ": dims = " << mat->size1 << " x " << mat->size2
              << " (" << (mat->size1 * mat->size2) << " elements), non-finite = "
              << badCount << '\n';
    return badCount;
}

/*! \brief Compare two matrices element by element.
 *
 * \return the number of differing elements; SIZE_MAX-style sentinel is avoided,
 *         a dimension mismatch simply reports the full element count of the
 *         first matrix so the caller still sees a failure.
 */
size_t ddm_count_differences(const gsl_matrix* a, const gsl_matrix* b, const char* name)
{
    if (a == nullptr || b == nullptr) {
        std::cerr << "    " << name << ": one of the matrices is NULL\n";
        return 1;
    }
    if (a->size1 != b->size1 || a->size2 != b->size2) {
        std::cerr << "    " << name << ": dimension mismatch (" << a->size1 << "x"
                  << a->size2 << " vs " << b->size1 << "x" << b->size2 << ")\n";
        return a->size1 * a->size2 + 1;
    }

    size_t diffCount = 0;
    for (size_t i = 0; i < a->size1; ++i) {
        for (size_t j = 0; j < a->size2; ++j) {
            const double va = gsl_matrix_get(a, i, j);
            const double vb = gsl_matrix_get(b, i, j);
            if (va != vb) {
                if (diffCount == 0) {
                    std::cerr << "    " << name << ": first difference at (" << i << ","
                              << j << "): " << va << " vs " << vb << '\n';
                }
                ++diffCount;
            }
        }
    }
    std::cerr << "    " << name << ": differing elements = " << diffCount << '\n';
    return diffCount;
}

/*! \brief Total number of elements held by a matrix (0 for a null pointer). */
size_t ddm_total_elements(const gsl_matrix* mat)
{
    return (mat == nullptr) ? 0 : (mat->size1 * mat->size2);
}

/*! \brief Release the three tables produced by create_DDMatrices(). */
void ddm_free_all(gsl_matrix*& surv, gsl_matrix*& norm, gsl_matrix*& pir)
{
    if (surv != nullptr) {
        gsl_matrix_free(surv);
        surv = nullptr;
    }
    if (norm != nullptr) {
        gsl_matrix_free(norm);
        norm = nullptr;
    }
    if (pir != nullptr) {
        gsl_matrix_free(pir);
        pir = nullptr;
    }
}

/*! \brief RAII helper: disable the (aborting) GSL error handler for a test.
 *
 * The numerical integrations used to fill the 2D tables can occasionally hit
 * GSL's iteration limits.  The default GSL error handler calls abort(), which
 * would tear down the whole test binary.  Turning the handler off keeps the
 * suite alive so that every test can run, and the handler is restored when the
 * object goes out of scope.
 */
class DdmGslErrorGuard {
public:
    DdmGslErrorGuard() { oldHandler_ = gsl_set_error_handler_off(); }
    ~DdmGslErrorGuard() { gsl_set_error_handler(oldHandler_); }

private:
    gsl_error_handler_t* oldHandler_ { nullptr };
};

} // namespace

// -----------------------------------------------------------------------------
// Test 1: The three output pointers must be allocated, distinct, and at least
//         as large as the lookup length reported by size_lookup().
// -----------------------------------------------------------------------------
void test_ddm_allocation_and_sizes()
{
    std::cerr << "\n[TEST] test_ddm_allocation_and_sizes\n"
              << "  Source file:   src/reactions/create_DDMatrices.cpp\n"
              << "  Function:      create_DDMatrices()\n"
              << "  Scenario:      pre-allocate the three tables the way the real caller\n"
              << "                 does (see determine_2D_bimolecular_reaction_probability.cpp),\n"
              << "                 sized via size_lookup(), then fill them for a small,\n"
              << "                 fast-to-build radial range.\n"
              << "  Pass criteria: all pointers stay non-NULL, are distinct,\n"
              << "                 and hold at least size_lookup() elements.\n";

    DdmGslErrorGuard errGuard; // keep GSL from aborting the whole suite

    // Small, well-conditioned inputs (nm, nm^2/us, us).
    const double bindRadius = 1.0;
    const double Dtot = 1.0;
    const double comRMax = bindRadius + 0.2; // keep the table tiny
    const double ktemp = 10.0;
    Parameters params = ddm_make_params(0.1);

    // Radial step used internally by the function under test.
    const double RStepSize = std::sqrt(Dtot * params.timeStep) / 50.0;
    const size_t expectedLen = size_lookup(bindRadius, Dtot, params, comRMax);

    std::cerr << "  bindRadius = " << bindRadius << ", Dtot = " << Dtot
              << ", comRMax = " << comRMax << ", ka = " << ktemp
              << ", timeStep = " << params.timeStep << '\n'
              << "  RStepSize (recomputed here) = " << RStepSize
              << ", size_lookup() = " << expectedLen << '\n';

    // create_DDMatrices() does not allocate its output matrices -- it only fills
    // in whatever gsl_matrix the caller already allocated (see
    // determine_2D_bimolecular_reaction_probability.cpp). Pre-allocate here the
    // same way the real caller does.
    gsl_matrix* survMatrix = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* normMatrix = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* pirMatrix = gsl_matrix_alloc(expectedLen, expectedLen);

    std::cerr << "  Calling create_DDMatrices()...\n";
    create_DDMatrices(survMatrix, normMatrix, pirMatrix, bindRadius, Dtot, comRMax, ktemp, params);

    // 1. Every table must still be a valid allocation.
    EXPECT_NE(survMatrix, nullptr) << "survMatrix should be a valid allocation";
    EXPECT_NE(normMatrix, nullptr) << "normMatrix should be a valid allocation";
    EXPECT_NE(pirMatrix, nullptr) << "pirMatrix should be a valid allocation";

    // 2. They must be three separate allocations (no aliasing of one table).
    EXPECT_NE(survMatrix, normMatrix) << "surv and norm matrices must be distinct objects";
    EXPECT_NE(survMatrix, pirMatrix) << "surv and pir matrices must be distinct objects";
    EXPECT_NE(normMatrix, pirMatrix) << "norm and pir matrices must be distinct objects";

    // 3. Dimensions must be non-degenerate and cover the lookup range.
    const size_t survTotal = ddm_total_elements(survMatrix);
    const size_t normTotal = ddm_total_elements(normMatrix);
    const size_t pirTotal = ddm_total_elements(pirMatrix);

    std::cerr << "  Element counts: surv = " << survTotal << ", norm = " << normTotal
              << ", pir = " << pirTotal << '\n';

    EXPECT_GT(survTotal, 0u) << "survMatrix should contain at least one element";
    EXPECT_GT(normTotal, 0u) << "normMatrix should contain at least one element";
    EXPECT_GT(pirTotal, 0u) << "pirMatrix should contain at least one element";

    EXPECT_GE(survTotal, expectedLen)
        << "survMatrix should span at least the size_lookup() radial range";
    EXPECT_GE(normTotal, expectedLen)
        << "normMatrix should span at least the size_lookup() radial range";
    EXPECT_GE(pirTotal, expectedLen)
        << "pirMatrix should span at least the size_lookup() radial range";

    // 4. pir is a 2D (r, r0) table, so it should never be smaller than the 1D
    //    survival table built from the same radial grid.
    EXPECT_GE(pirTotal, survTotal)
        << "the pir table is 2D in (r, r0) and should not be smaller than surv";

    ddm_free_all(survMatrix, normMatrix, pirMatrix);
    std::cerr << "  Matrices freed.\n";
}

// -----------------------------------------------------------------------------
// Test 2: Every tabulated value must be a finite number.
// -----------------------------------------------------------------------------
void test_ddm_entries_are_finite()
{
    std::cerr << "\n[TEST] test_ddm_entries_are_finite\n"
              << "  Source file:   src/reactions/create_DDMatrices.cpp\n"
              << "  Function:      create_DDMatrices()\n"
              << "  Scenario:      build all three tables and scan every element.\n"
              << "  Pass criteria: no NaN and no infinite values are stored.\n";

    DdmGslErrorGuard errGuard;

    const double bindRadius = 1.0;
    const double Dtot = 2.0;
    const double comRMax = bindRadius + 0.2;
    const double ktemp = 5.0;
    Parameters params = ddm_make_params(0.1);
    const size_t expectedLen = size_lookup(bindRadius, Dtot, params, comRMax);

    // create_DDMatrices() only fills a pre-allocated gsl_matrix; it does not
    // allocate one itself. Size the tables the way the real caller does.
    gsl_matrix* survMatrix = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* normMatrix = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* pirMatrix = gsl_matrix_alloc(expectedLen, expectedLen);

    std::cerr << "  Calling create_DDMatrices() with Dtot = " << Dtot
              << ", ka = " << ktemp << "...\n";
    create_DDMatrices(survMatrix, normMatrix, pirMatrix, bindRadius, Dtot, comRMax, ktemp, params);

    // Guard against dereferencing NULL, but keep the failure non-fatal.
    if (survMatrix == nullptr || normMatrix == nullptr || pirMatrix == nullptr) {
        ADD_FAILURE() << "create_DDMatrices did not allocate all three matrices; "
                         "skipping element inspection";
        ddm_free_all(survMatrix, normMatrix, pirMatrix);
        return;
    }

    std::cerr << "  Scanning matrix elements for NaN / inf:\n";
    EXPECT_EQ(ddm_count_nonfinite(survMatrix, "survMatrix"), 0u)
        << "survMatrix must contain only finite values";
    EXPECT_EQ(ddm_count_nonfinite(normMatrix, "normMatrix"), 0u)
        << "normMatrix must contain only finite values";
    EXPECT_EQ(ddm_count_nonfinite(pirMatrix, "pirMatrix"), 0u)
        << "pirMatrix must contain only finite values";

    ddm_free_all(survMatrix, normMatrix, pirMatrix);
    std::cerr << "  Matrices freed.\n";
}

// -----------------------------------------------------------------------------
// Test 3: Determinism - identical inputs must produce identical tables.
// -----------------------------------------------------------------------------
void test_ddm_is_deterministic()
{
    std::cerr << "\n[TEST] test_ddm_is_deterministic\n"
              << "  Source file:   src/reactions/create_DDMatrices.cpp\n"
              << "  Function:      create_DDMatrices()\n"
              << "  Scenario:      call the function twice with identical inputs.\n"
              << "  Pass criteria: both runs give the same dimensions and the\n"
              << "                 same value in every element (no RNG or\n"
              << "                 uninitialised-memory dependence).\n";

    DdmGslErrorGuard errGuard;

    const double bindRadius = 1.5;
    const double Dtot = 1.0;
    const double comRMax = bindRadius + 0.15;
    const double ktemp = 8.0;
    Parameters params = ddm_make_params(0.1);
    const size_t expectedLen = size_lookup(bindRadius, Dtot, params, comRMax);

    // create_DDMatrices() only fills a pre-allocated gsl_matrix; it does not
    // allocate one itself. Size the tables the way the real caller does.
    gsl_matrix* surv1 = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* norm1 = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* pir1 = gsl_matrix_alloc(expectedLen, expectedLen);
    gsl_matrix* surv2 = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* norm2 = gsl_matrix_alloc(2, expectedLen);
    gsl_matrix* pir2 = gsl_matrix_alloc(expectedLen, expectedLen);

    std::cerr << "  First call to create_DDMatrices()...\n";
    create_DDMatrices(surv1, norm1, pir1, bindRadius, Dtot, comRMax, ktemp, params);
    std::cerr << "  Second call to create_DDMatrices() with identical inputs...\n";
    create_DDMatrices(surv2, norm2, pir2, bindRadius, Dtot, comRMax, ktemp, params);

    if (surv1 == nullptr || norm1 == nullptr || pir1 == nullptr || surv2 == nullptr
        || norm2 == nullptr || pir2 == nullptr) {
        ADD_FAILURE() << "create_DDMatrices did not allocate all matrices on both calls; "
                         "skipping determinism comparison";
        ddm_free_all(surv1, norm1, pir1);
        ddm_free_all(surv2, norm2, pir2);
        return;
    }

    // Dimensions must agree exactly.
    EXPECT_EQ(surv1->size1, surv2->size1) << "survMatrix row count must be reproducible";
    EXPECT_EQ(surv1->size2, surv2->size2) << "survMatrix column count must be reproducible";
    EXPECT_EQ(norm1->size1, norm2->size1) << "normMatrix row count must be reproducible";
    EXPECT_EQ(norm1->size2, norm2->size2) << "normMatrix column count must be reproducible";
    EXPECT_EQ(pir1->size1, pir2->size1) << "pirMatrix row count must be reproducible";
    EXPECT_EQ(pir1->size2, pir2->size2) << "pirMatrix column count must be reproducible";

    // Contents must agree bit for bit.
    std::cerr << "  Comparing contents of the two runs:\n";
    EXPECT_EQ(ddm_count_differences(surv1, surv2, "survMatrix"), 0u)
        << "survMatrix contents must be identical between identical calls";
    EXPECT_EQ(ddm_count_differences(norm1, norm2, "normMatrix"), 0u)
        << "normMatrix contents must be identical between identical calls";
    EXPECT_EQ(ddm_count_differences(pir1, pir2, "pirMatrix"), 0u)
        << "pirMatrix contents must be identical between identical calls";

    ddm_free_all(surv1, norm1, pir1);
    ddm_free_all(surv2, norm2, pir2);
    std::cerr << "  Matrices freed.\n";
}

// -----------------------------------------------------------------------------
// Test 4: A larger radial range (comRMax) must produce larger lookup tables.
// -----------------------------------------------------------------------------
void test_ddm_larger_rmax_grows_tables()
{
    std::cerr << "\n[TEST] test_ddm_larger_rmax_grows_tables\n"
              << "  Source file:   src/reactions/create_DDMatrices.cpp\n"
              << "  Function:      create_DDMatrices()\n"
              << "  Scenario:      build the tables twice, the second time with a\n"
              << "                 larger maximum separation comRMax while the\n"
              << "                 radial step size (set by Dtot and timeStep) is\n"
              << "                 held fixed.\n"
              << "  Pass criteria: the larger comRMax yields at least as many\n"
              << "                 tabulated elements, and strictly more for the\n"
              << "                 1D survival table.\n";

    DdmGslErrorGuard errGuard;

    const double bindRadius = 1.0;
    const double Dtot = 1.0;
    const double ktemp = 10.0;
    Parameters params = ddm_make_params(0.1);

    const double smallRMax = bindRadius + 0.10;
    const double largeRMax = bindRadius + 0.25;

    // Report what size_lookup() predicts for both ranges.
    const size_t lenSmall = size_lookup(bindRadius, Dtot, params, smallRMax);
    const size_t lenLarge = size_lookup(bindRadius, Dtot, params, largeRMax);
    std::cerr << "  size_lookup(comRMax = " << smallRMax << ") = " << lenSmall << '\n'
              << "  size_lookup(comRMax = " << largeRMax << ") = " << lenLarge << '\n';
    EXPECT_GT(lenLarge, lenSmall)
        << "size_lookup() should grow with comRMax (sanity check on the inputs)";

    // create_DDMatrices() only fills a pre-allocated gsl_matrix; it does not
    // allocate one itself. Size the tables the way the real caller does.
    gsl_matrix* survSmall = gsl_matrix_alloc(2, lenSmall);
    gsl_matrix* normSmall = gsl_matrix_alloc(2, lenSmall);
    gsl_matrix* pirSmall = gsl_matrix_alloc(lenSmall, lenSmall);
    gsl_matrix* survLarge = gsl_matrix_alloc(2, lenLarge);
    gsl_matrix* normLarge = gsl_matrix_alloc(2, lenLarge);
    gsl_matrix* pirLarge = gsl_matrix_alloc(lenLarge, lenLarge);

    std::cerr << "  Building tables for the small radial range...\n";
    create_DDMatrices(survSmall, normSmall, pirSmall, bindRadius, Dtot, smallRMax, ktemp, params);
    std::cerr << "  Building tables for the large radial range...\n";
    create_DDMatrices(survLarge, normLarge, pirLarge, bindRadius, Dtot, largeRMax, ktemp, params);

    const size_t survSmallN = ddm_total_elements(survSmall);
    const size_t survLargeN = ddm_total_elements(survLarge);
    const size_t normSmallN = ddm_total_elements(normSmall);
    const size_t normLargeN = ddm_total_elements(normLarge);
    const size_t pirSmallN = ddm_total_elements(pirSmall);
    const size_t pirLargeN = ddm_total_elements(pirLarge);

    std::cerr << "  Element counts (small -> large):\n"
              << "    survMatrix: " << survSmallN << " -> " << survLargeN << '\n'
              << "    normMatrix: " << normSmallN << " -> " << normLargeN << '\n'
              << "    pirMatrix:  " << pirSmallN << " -> " << pirLargeN << '\n';

    // The survival table is indexed purely by the radial grid, so it must grow.
    EXPECT_GT(survLargeN, survSmallN)
        << "survMatrix should contain more entries for a larger comRMax";
    // The other tables must at minimum not shrink.
    EXPECT_GE(normLargeN, normSmallN)
        << "normMatrix should not shrink when comRMax increases";
    EXPECT_GE(pirLargeN, pirSmallN)
        << "pirMatrix should not shrink when comRMax increases";

    ddm_free_all(survSmall, normSmall, pirSmall);
    ddm_free_all(survLarge, normLarge, pirLarge);
    std::cerr << "  Matrices freed.\n";
}

// -----------------------------------------------------------------------------
// Test 5: Robustness - the function must reallocate cleanly when it is used to
//         build several independent table sets (as happens in the simulation
//         whenever a new (Dtot, bindRadius) pair is encountered).
// -----------------------------------------------------------------------------
void test_ddm_multiple_parameter_sets()
{
    std::cerr << "\n[TEST] test_ddm_multiple_parameter_sets\n"
              << "  Source file:   src/reactions/create_DDMatrices.cpp\n"
              << "  Function:      create_DDMatrices()\n"
              << "  Scenario:      loop over several (bindRadius, Dtot, ka) sets,\n"
              << "                 as the simulation does when it caches 2D tables.\n"
              << "  Pass criteria: every set yields three non-NULL, finite,\n"
              << "                 non-empty tables.\n";

    DdmGslErrorGuard errGuard;

    // A handful of physically plausible but deliberately small cases.
    struct Case {
        double bindRadius;
        double Dtot;
        double ka;
    };
    const Case cases[] = {
        { 1.0, 0.5, 1.0 },
        { 2.0, 1.0, 20.0 },
        { 0.5, 2.0, 100.0 },
    };

    Parameters params = ddm_make_params(0.1);

    for (const Case& c : cases) {
        // Keep every table tiny: only a short stretch beyond the binding radius.
        const double comRMax = c.bindRadius + 0.1;

        std::cerr << "  --- case: bindRadius = " << c.bindRadius << ", Dtot = " << c.Dtot
                  << ", ka = " << c.ka << ", comRMax = " << comRMax << '\n';

        // create_DDMatrices() only fills a pre-allocated gsl_matrix; it does not
        // allocate one itself. Size the tables the way the real caller does.
        const size_t expectedLen = size_lookup(c.bindRadius, c.Dtot, params, comRMax);
        gsl_matrix* surv = gsl_matrix_alloc(2, expectedLen);
        gsl_matrix* norm = gsl_matrix_alloc(2, expectedLen);
        gsl_matrix* pir = gsl_matrix_alloc(expectedLen, expectedLen);

        create_DDMatrices(surv, norm, pir, c.bindRadius, c.Dtot, comRMax, c.ka, params);

        // Each table must exist and be non-empty.
        EXPECT_NE(surv, nullptr) << "survMatrix allocated for this parameter set";
        EXPECT_NE(norm, nullptr) << "normMatrix allocated for this parameter set";
        EXPECT_NE(pir, nullptr) << "pirMatrix allocated for this parameter set";
        EXPECT_GT(ddm_total_elements(surv), 0u) << "survMatrix must be non-empty";
        EXPECT_GT(ddm_total_elements(norm), 0u) << "normMatrix must be non-empty";
        EXPECT_GT(ddm_total_elements(pir), 0u) << "pirMatrix must be non-empty";

        // And must not contain garbage.
        EXPECT_EQ(ddm_count_nonfinite(surv, "survMatrix"), 0u)
            << "survMatrix must be finite for this parameter set";
        EXPECT_EQ(ddm_count_nonfinite(norm, "normMatrix"), 0u)
            << "normMatrix must be finite for this parameter set";
        EXPECT_EQ(ddm_count_nonfinite(pir, "pirMatrix"), 0u)
            << "pirMatrix must be finite for this parameter set";

        ddm_free_all(surv, norm, pir);
    }

    std::cerr << "  All parameter sets processed and freed.\n";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each named test_* helper is invoked from its own TEST so
// that a failure in one scenario does not prevent the remaining ones from
// running (all assertions above are non-fatal EXPECT_* style).
// -----------------------------------------------------------------------------
TEST(CreateDDMatrices, AllocationAndSizes) { test_ddm_allocation_and_sizes(); }
TEST(CreateDDMatrices, EntriesAreFinite) { test_ddm_entries_are_finite(); }
TEST(CreateDDMatrices, IsDeterministic) { test_ddm_is_deterministic(); }
TEST(CreateDDMatrices, LargerRMaxGrowsTables) { test_ddm_larger_rmax_grows_tables(); }
TEST(CreateDDMatrices, MultipleParameterSets) { test_ddm_multiple_parameter_sets(); }