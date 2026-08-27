/*! \file test_rand_gsl.cpp
 *
 * ### Unit tests for ../src/math/rand_gsl.cpp
 *
 * The file under test is a thin wrapper around the GSL random number
 * generator.  It provides:
 *
 *   - double rand_gsl()                  : uniform double in [0,1)
 *   - double rand_gsl64()                : uniform double with extra low bits
 *   - void   srand_gsl(int)              : re-seeds the (already allocated) RNG
 *   - void   write_rng_state()           : dumps RNG state to "DATA/rng_state"
 *   - void   write_rng_state(int rank)   : dumps to "DATA/rng_state_<rank>"
 *   - void   write_rng_state_simItr(int) : dumps to "RESTARTS/rng_state<itr>"
 *   - void   read_rng_state()            : restores from "./rng_state"
 *   - void   read_rng_state(int rank)    : restores from "./rng_state_<rank>"
 *   - double GaussV()                    : Box-Muller (polar) normal deviate
 *
 * Notes on things deliberately NOT exercised here (they would kill the whole
 * gtest binary rather than fail a single case):
 *
 *   - write_rng_state*() call ferror() on the FILE* returned by fopen() without
 *     first checking for nullptr.  If the target directory does not exist the
 *     pointer is null and ferror(nullptr) is undefined behaviour.  Every test
 *     below therefore creates "DATA/" and "RESTARTS/" first, and skips the
 *     assertions (with a warning) if the directory cannot be created.
 *   - read_rng_state(int rank) calls fclose(stateIn) on the *failure* path,
 *     where stateIn may be nullptr.  We therefore only call the rank overload
 *     with a file that is guaranteed to exist.  (The no-argument
 *     read_rng_state() returns before fclose() on failure, so the missing-file
 *     branch of that overload *is* safe and is tested.)
 *   - srand_gsl() is not called: the suite convention is to allocate/seed the
 *     generator directly with gsl_rng_alloc()/gsl_rng_set().  srand_gsl(num) is
 *     literally "gsl_rng_set(r, num)", and the equivalent behaviour (identical
 *     streams for identical seeds, different streams for different seeds) is
 *     verified through gsl_rng_set() below.
 */

#include "math/rand_gsl.hpp"

#include <gtest/gtest.h>

#include <sys/stat.h>

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

// The global generator is defined once, in unit_tests/src/gtest_main.cpp.
extern gsl_rng* r;

namespace {

// -----------------------------------------------------------------------------
// Helper: make sure the global generator exists and is seeded reproducibly.
//
// Per the suite convention we never call srand_gsl() for initialisation, since
// it only re-seeds and would dereference a null generator pointer.  We allocate
// with the default GSL type and seed with 42 so that every test starts from a
// known, repeatable stream.
// -----------------------------------------------------------------------------
void randgsl_init_rng(unsigned long seed = 42)
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
    }
    gsl_rng_set(r, seed);
}

/*! \brief Create a directory if it does not already exist.
 *  \return true if the directory exists (or was created) afterwards.
 */
bool randgsl_ensure_dir(const char* dirName)
{
    struct stat info;
    if (stat(dirName, &info) == 0)
        return true; // already there (file or directory)
    mkdir(dirName, 0755);
    return stat(dirName, &info) == 0;
}

/*! \brief Return the size of a file in bytes, or -1 if it cannot be opened. */
long randgsl_file_size(const char* path)
{
    FILE* f = fopen(path, "rb");
    if (f == nullptr)
        return -1;
    fseek(f, 0, SEEK_END);
    long sz = ftell(f);
    fclose(f);
    return sz;
}

/*! \brief Byte-for-byte copy of a file. \return true on success. */
bool randgsl_copy_file(const char* src, const char* dst)
{
    FILE* in = fopen(src, "rb");
    if (in == nullptr)
        return false;
    FILE* out = fopen(dst, "wb");
    if (out == nullptr) {
        fclose(in);
        return false;
    }
    char buf[4096];
    size_t n = 0;
    while ((n = fread(buf, 1, sizeof(buf), in)) > 0)
        fwrite(buf, 1, n, out);
    fclose(in);
    fclose(out);
    return true;
}

} // namespace

// -----------------------------------------------------------------------------
// rand_gsl(): range, determinism for a fixed seed, and variability.
// -----------------------------------------------------------------------------
void test_randgsl_rand_gsl_basic()
{
    std::cerr << "\n[TEST] test_randgsl_rand_gsl_basic\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      rand_gsl()\n"
              << "  Checks:        (a) every draw lies in [0,1)\n"
              << "                 (b) the same seed reproduces the same stream\n"
              << "                 (c) a different seed gives a different stream\n"
              << "                 (d) draws are not all identical\n";

    const int kNumDraws = 5000;

    // (a) + (d): range check and simple "is it actually varying" check.
    randgsl_init_rng(42);
    double firstDraw = rand_gsl();
    bool allInRange = true;
    bool sawDifferentValue = false;
    for (int i = 1; i < kNumDraws; ++i) {
        double v = rand_gsl();
        if (!(v >= 0.0 && v < 1.0))
            allInRange = false;
        if (v != firstDraw)
            sawDifferentValue = true;
    }
    std::cerr << "  -> Drew " << kNumDraws << " values, first = " << firstDraw << "\n";
    EXPECT_TRUE(allInRange) << "rand_gsl() must always return a value in [0,1)";
    EXPECT_TRUE(sawDifferentValue) << "rand_gsl() returned the same value every time";

    // (b) determinism: re-seeding with 42 must replay the identical sequence.
    randgsl_init_rng(42);
    std::vector<double> seqA;
    for (int i = 0; i < 10; ++i)
        seqA.push_back(rand_gsl());

    randgsl_init_rng(42);
    std::vector<double> seqB;
    for (int i = 0; i < 10; ++i)
        seqB.push_back(rand_gsl());

    std::cerr << "  -> Comparing two 10-value streams generated from seed 42\n";
    for (size_t i = 0; i < seqA.size(); ++i) {
        EXPECT_DOUBLE_EQ(seqA[i], seqB[i])
            << "Draw " << i << " differs between two runs seeded identically";
    }

    // (c) a different seed should (with overwhelming probability) differ.
    randgsl_init_rng(12345);
    std::vector<double> seqC;
    for (int i = 0; i < 10; ++i)
        seqC.push_back(rand_gsl());

    bool anyDifferent = false;
    for (size_t i = 0; i < seqA.size(); ++i) {
        if (seqA[i] != seqC[i])
            anyDifferent = true;
    }
    std::cerr << "  -> Comparing seed 42 stream against seed 12345 stream\n";
    EXPECT_TRUE(anyDifferent) << "Different seeds produced an identical 10-value stream";
}

// -----------------------------------------------------------------------------
// rand_gsl(): coarse uniformity check (mean and bin occupancy).
// -----------------------------------------------------------------------------
void test_randgsl_rand_gsl_uniformity()
{
    std::cerr << "\n[TEST] test_randgsl_rand_gsl_uniformity\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      rand_gsl()\n"
              << "  Checks:        sample mean ~= 0.5 and all 10 equal-width bins\n"
              << "                 are populated (loose statistical sanity test).\n";

    const int kNumDraws = 100000;
    randgsl_init_rng(42);

    double sum = 0.0;
    std::vector<int> bins(10, 0);
    for (int i = 0; i < kNumDraws; ++i) {
        double v = rand_gsl();
        sum += v;
        int b = static_cast<int>(v * 10.0);
        if (b < 0)
            b = 0;
        if (b > 9)
            b = 9;
        ++bins[b];
    }
    double mean = sum / kNumDraws;
    std::cerr << "  -> Sample mean over " << kNumDraws << " draws = " << mean
              << " (expected 0.5)\n";

    // Standard error of the mean for U(0,1) with 1e5 samples is ~9.1e-4, so a
    // 0.01 window is ~11 sigma: it will not fire by chance.
    EXPECT_NEAR(mean, 0.5, 0.01) << "rand_gsl() mean is far from the expected 0.5";

    // Each bin should hold roughly 10000 counts; allow a very generous window.
    for (size_t b = 0; b < bins.size(); ++b) {
        std::cerr << "     bin[" << b << "] = " << bins[b] << "\n";
        EXPECT_GT(bins[b], kNumDraws / 20)
            << "Bin " << b << " is badly under-populated; distribution is not uniform";
        EXPECT_LT(bins[b], kNumDraws / 5)
            << "Bin " << b << " is badly over-populated; distribution is not uniform";
    }
}

// -----------------------------------------------------------------------------
// rand_gsl64(): range, determinism, and the fact that it consumes two draws.
// -----------------------------------------------------------------------------
void test_randgsl_rand_gsl64()
{
    std::cerr << "\n[TEST] test_randgsl_rand_gsl64\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      rand_gsl64()\n"
              << "  Implementation: u1 + u2/(gsl_rng_max+1), i.e. a low-order\n"
              << "                  refinement of a uniform draw.\n"
              << "  Checks:        (a) result >= 0 and < 1 + 1/(max+1)\n"
              << "                 (b) identical seeds reproduce the stream\n"
              << "                 (c) it consumes exactly two underlying draws\n";

    randgsl_init_rng(42);
    const double tinyBit = 1.0 / (static_cast<double>(gsl_rng_max(r)) + 1.0);
    std::cerr << "  -> 1/(gsl_rng_max+1) = " << tinyBit << "\n";

    // (a) range.
    bool allInRange = true;
    for (int i = 0; i < 5000; ++i) {
        double v = rand_gsl64();
        if (!(v >= 0.0 && v < 1.0 + tinyBit))
            allInRange = false;
    }
    EXPECT_TRUE(allInRange)
        << "rand_gsl64() returned a value outside [0, 1 + 1/(max+1))";

    // (b) determinism.
    randgsl_init_rng(7);
    std::vector<double> seqA;
    for (int i = 0; i < 8; ++i)
        seqA.push_back(rand_gsl64());
    randgsl_init_rng(7);
    std::vector<double> seqB;
    for (int i = 0; i < 8; ++i)
        seqB.push_back(rand_gsl64());
    std::cerr << "  -> Comparing two 8-value rand_gsl64 streams from seed 7\n";
    for (size_t i = 0; i < seqA.size(); ++i)
        EXPECT_DOUBLE_EQ(seqA[i], seqB[i]) << "rand_gsl64 draw " << i << " is not reproducible";

    // (c) rand_gsl64 is u1 + tinyBit*u2, so replaying the raw generator by hand
    //     from the same seed must reproduce it exactly (proving two draws are used).
    randgsl_init_rng(99);
    double fromWrapper = rand_gsl64();
    randgsl_init_rng(99);
    double u1 = gsl_rng_uniform(r);
    double u2 = gsl_rng_uniform(r);
    double byHand = u1 + tinyBit * u2;
    std::cerr << "  -> wrapper = " << fromWrapper << ", hand-computed u1 + u2/(max+1) = "
              << byHand << "\n";
    EXPECT_DOUBLE_EQ(fromWrapper, byHand)
        << "rand_gsl64() does not equal u1 + u2/(gsl_rng_max+1) (draw count mismatch?)";
}

// -----------------------------------------------------------------------------
// srand_gsl(): documented, and the equivalent gsl_rng_set() behaviour verified.
// -----------------------------------------------------------------------------
void test_randgsl_seeding_semantics()
{
    std::cerr << "\n[TEST] test_randgsl_seeding_semantics\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      srand_gsl(int)  [NOT called - see note]\n"
              << "  Note:          srand_gsl(num) is exactly gsl_rng_set(r, num);\n"
              << "                 it never allocates the generator, so the suite\n"
              << "                 initialises with gsl_rng_alloc + gsl_rng_set.\n"
              << "  Checks:        re-seeding rewinds the stream, and two distinct\n"
              << "                 seeds diverge.\n";

    // Seeding with the same value must rewind the stream to the same point.
    randgsl_init_rng(2024);
    double a1 = rand_gsl();
    double a2 = rand_gsl();
    randgsl_init_rng(2024);
    double b1 = rand_gsl();
    double b2 = rand_gsl();
    std::cerr << "  -> seed 2024 draws: (" << a1 << ", " << a2 << ") vs ("
              << b1 << ", " << b2 << ")\n";
    EXPECT_DOUBLE_EQ(a1, b1) << "Re-seeding did not rewind the generator";
    EXPECT_DOUBLE_EQ(a2, b2) << "Re-seeding did not rewind the generator";

    // Distinct seeds must diverge.
    randgsl_init_rng(1);
    double c1 = rand_gsl();
    randgsl_init_rng(2);
    double d1 = rand_gsl();
    std::cerr << "  -> seed 1 first draw = " << c1 << ", seed 2 first draw = " << d1 << "\n";
    EXPECT_NE(c1, d1) << "Two different seeds produced the same first draw";
}

// -----------------------------------------------------------------------------
// GaussV(): finiteness, determinism, and normal-distribution statistics.
// -----------------------------------------------------------------------------
void test_randgsl_gaussv()
{
    std::cerr << "\n[TEST] test_randgsl_gaussv\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      GaussV()\n"
              << "  Checks:        (a) all returns are finite\n"
              << "                 (b) the same seed reproduces the same stream\n"
              << "                 (c) sample mean ~ 0 and sample sd ~ 1\n"
              << "                 (d) roughly 2/3 of samples fall within +/-1 sd\n";

    // (a)+(c)+(d): statistics over a large sample.
    const int kNumSamples = 50000;
    randgsl_init_rng(42);

    double sum = 0.0;
    double sumSq = 0.0;
    int withinOneSigma = 0;
    bool allFinite = true;
    for (int i = 0; i < kNumSamples; ++i) {
        double g = GaussV();
        if (!std::isfinite(g))
            allFinite = false;
        sum += g;
        sumSq += g * g;
        if (std::fabs(g) <= 1.0)
            ++withinOneSigma;
    }
    double mean = sum / kNumSamples;
    double variance = (sumSq / kNumSamples) - (mean * mean);
    double sd = std::sqrt(variance);
    double fracWithin = static_cast<double>(withinOneSigma) / kNumSamples;

    std::cerr << "  -> " << kNumSamples << " samples: mean = " << mean
              << ", sd = " << sd << ", fraction within +/-1 = " << fracWithin << "\n";

    EXPECT_TRUE(allFinite) << "GaussV() returned a non-finite value";
    // SEM for N(0,1) with 5e4 samples is ~0.0045, so 0.05 is >10 sigma.
    EXPECT_NEAR(mean, 0.0, 0.05) << "GaussV() mean is not consistent with 0";
    EXPECT_NEAR(sd, 1.0, 0.05) << "GaussV() standard deviation is not consistent with 1";
    // Exact value for a standard normal is 0.6827.
    EXPECT_NEAR(fracWithin, 0.6827, 0.02)
        << "Fraction of GaussV() samples within +/-1 sd is not normal-like";

    // (b) determinism: identical seeds must reproduce the deviates exactly.
    randgsl_init_rng(2718);
    std::vector<double> seqA;
    for (int i = 0; i < 6; ++i)
        seqA.push_back(GaussV());
    randgsl_init_rng(2718);
    std::vector<double> seqB;
    for (int i = 0; i < 6; ++i)
        seqB.push_back(GaussV());

    std::cerr << "  -> Comparing two 6-value GaussV streams from seed 2718\n";
    for (size_t i = 0; i < seqA.size(); ++i) {
        std::cerr << "     [" << i << "] " << seqA[i] << " vs " << seqB[i] << "\n";
        EXPECT_DOUBLE_EQ(seqA[i], seqB[i]) << "GaussV draw " << i << " is not reproducible";
    }
}

// -----------------------------------------------------------------------------
// write_rng_state(): writes "DATA/rng_state" with gsl_rng_size(r) bytes.
// -----------------------------------------------------------------------------
void test_randgsl_write_rng_state()
{
    std::cerr << "\n[TEST] test_randgsl_write_rng_state\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      write_rng_state()\n"
              << "  Checks:        the file DATA/rng_state is created and its size\n"
              << "                 equals gsl_rng_size(r).\n"
              << "  Note:          DATA/ is created first - the function calls\n"
              << "                 ferror() on a possibly-null FILE*.\n";

    randgsl_init_rng(42);

    if (!randgsl_ensure_dir("DATA")) {
        std::cerr << "  !! Could not create DATA/ - skipping write_rng_state() assertions\n";
        ADD_FAILURE() << "Could not create the DATA directory needed by write_rng_state()";
        return;
    }

    remove("DATA/rng_state"); // start clean so the check below is meaningful
    write_rng_state();

    long sz = randgsl_file_size("DATA/rng_state");
    long expected = static_cast<long>(gsl_rng_size(r));
    std::cerr << "  -> DATA/rng_state size = " << sz << " bytes (gsl_rng_size = "
              << expected << ")\n";

    EXPECT_GT(sz, 0) << "write_rng_state() did not create a non-empty DATA/rng_state";
    EXPECT_EQ(sz, expected)
        << "write_rng_state() wrote a file whose size differs from gsl_rng_size(r)";
}

// -----------------------------------------------------------------------------
// write_rng_state(int rank): writes "DATA/rng_state_<rank>".
// -----------------------------------------------------------------------------
void test_randgsl_write_rng_state_rank()
{
    std::cerr << "\n[TEST] test_randgsl_write_rng_state_rank\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      write_rng_state(int rank)\n"
              << "  Checks:        DATA/rng_state_<rank> is created, is non-empty,\n"
              << "                 and two different ranks produce separate files.\n";

    randgsl_init_rng(42);

    if (!randgsl_ensure_dir("DATA")) {
        std::cerr << "  !! Could not create DATA/ - skipping assertions\n";
        ADD_FAILURE() << "Could not create the DATA directory needed by write_rng_state(rank)";
        return;
    }

    const int rankA = 3;
    const int rankB = 7;
    const std::string fileA = "DATA/rng_state_" + std::to_string(rankA);
    const std::string fileB = "DATA/rng_state_" + std::to_string(rankB);

    remove(fileA.c_str());
    remove(fileB.c_str());

    write_rng_state(rankA);
    write_rng_state(rankB);

    long szA = randgsl_file_size(fileA.c_str());
    long szB = randgsl_file_size(fileB.c_str());
    long expected = static_cast<long>(gsl_rng_size(r));
    std::cerr << "  -> " << fileA << " size = " << szA << "\n"
              << "  -> " << fileB << " size = " << szB << " (expected " << expected << ")\n";

    EXPECT_EQ(szA, expected) << "write_rng_state(" << rankA << ") wrote an unexpected size";
    EXPECT_EQ(szB, expected) << "write_rng_state(" << rankB << ") wrote an unexpected size";
}

// -----------------------------------------------------------------------------
// write_rng_state_simItr(int): writes "RESTARTS/rng_state<simItr>".
// -----------------------------------------------------------------------------
void test_randgsl_write_rng_state_simItr()
{
    std::cerr << "\n[TEST] test_randgsl_write_rng_state_simItr\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      write_rng_state_simItr(int simItr)\n"
              << "  Checks:        RESTARTS/rng_state<simItr> is created with\n"
              << "                 gsl_rng_size(r) bytes.\n"
              << "  Note:          RESTARTS/ is created first (null-FILE* hazard).\n";

    randgsl_init_rng(42);

    if (!randgsl_ensure_dir("RESTARTS")) {
        std::cerr << "  !! Could not create RESTARTS/ - skipping assertions\n";
        ADD_FAILURE() << "Could not create the RESTARTS directory needed by write_rng_state_simItr()";
        return;
    }

    const int simItr = 17;
    const std::string fname = "RESTARTS/rng_state" + std::to_string(simItr);
    remove(fname.c_str());

    write_rng_state_simItr(simItr);

    long sz = randgsl_file_size(fname.c_str());
    long expected = static_cast<long>(gsl_rng_size(r));
    std::cerr << "  -> " << fname << " size = " << sz << " (expected " << expected << ")\n";

    EXPECT_GT(sz, 0) << "write_rng_state_simItr() did not create a non-empty file";
    EXPECT_EQ(sz, expected)
        << "write_rng_state_simItr() wrote a file whose size differs from gsl_rng_size(r)";
}

// -----------------------------------------------------------------------------
// read_rng_state(): round trip - saved state must reproduce the future stream.
// -----------------------------------------------------------------------------
void test_randgsl_read_rng_state_roundtrip()
{
    std::cerr << "\n[TEST] test_randgsl_read_rng_state_roundtrip\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Functions:     write_rng_state() + read_rng_state()\n"
              << "  Scenario:      save the state mid-stream, keep drawing, restore,\n"
              << "                 draw again.\n"
              << "  Pass criteria: the draws after the restore equal the draws that\n"
              << "                 immediately followed the save.\n";

    if (!randgsl_ensure_dir("DATA")) {
        std::cerr << "  !! Could not create DATA/ - skipping assertions\n";
        ADD_FAILURE() << "Could not create the DATA directory needed by write_rng_state()";
        return;
    }

    randgsl_init_rng(42);

    // Burn a few draws so the saved state is not simply "freshly seeded".
    for (int i = 0; i < 5; ++i)
        rand_gsl();

    // write_rng_state() always writes to DATA/rng_state, but read_rng_state()
    // reads from ./rng_state, so we copy the file into place.
    write_rng_state();
    ASSERT_TRUE(randgsl_copy_file("DATA/rng_state", "rng_state"))
        << "Could not stage DATA/rng_state as ./rng_state for read_rng_state()";

    // Record the "expected future" of the saved state.
    std::vector<double> expected;
    for (int i = 0; i < 5; ++i)
        expected.push_back(rand_gsl());

    // Advance the generator further so a no-op restore would be detected.
    for (int i = 0; i < 20; ++i)
        rand_gsl();

    std::cerr << "  -> Calling read_rng_state() to rewind to the saved state\n";
    read_rng_state();

    std::vector<double> restored;
    for (int i = 0; i < 5; ++i)
        restored.push_back(rand_gsl());

    for (size_t i = 0; i < expected.size(); ++i) {
        std::cerr << "     [" << i << "] expected " << expected[i]
                  << " vs restored " << restored[i] << "\n";
        EXPECT_DOUBLE_EQ(restored[i], expected[i])
            << "read_rng_state() did not restore the saved generator state (draw " << i << ")";
    }

    // Clean up the staged file so the "missing file" test below is deterministic.
    remove("rng_state");
}

// -----------------------------------------------------------------------------
// read_rng_state(): missing file must be handled gracefully (no state change).
// -----------------------------------------------------------------------------
void test_randgsl_read_rng_state_missing_file()
{
    std::cerr << "\n[TEST] test_randgsl_read_rng_state_missing_file\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Function:      read_rng_state()\n"
              << "  Scenario:      ./rng_state does not exist.\n"
              << "  Pass criteria: the call returns without touching the generator,\n"
              << "                 so the next draw matches the un-disturbed stream.\n"
              << "  Note:          only the no-argument overload is safe here; the\n"
              << "                 rank overload fcloses a null FILE* on failure.\n";

    // Make sure the file really is absent.
    remove("rng_state");

    // Reference stream from a fresh seed.
    randgsl_init_rng(4242);
    double reference = rand_gsl();

    // Same seed, but attempt a (failing) restore before drawing.
    randgsl_init_rng(4242);
    std::cerr << "  -> Calling read_rng_state() with no ./rng_state present\n";
    read_rng_state();
    double afterFailedRead = rand_gsl();

    std::cerr << "  -> reference draw = " << reference
              << ", draw after failed read = " << afterFailedRead << "\n";
    EXPECT_DOUBLE_EQ(afterFailedRead, reference)
        << "A failed read_rng_state() must leave the generator state untouched";
}

// -----------------------------------------------------------------------------
// read_rng_state(int rank): round trip through DATA/rng_state_<rank>.
// -----------------------------------------------------------------------------
void test_randgsl_read_rng_state_rank_roundtrip()
{
    std::cerr << "\n[TEST] test_randgsl_read_rng_state_rank_roundtrip\n"
              << "  Source file:   src/math/rand_gsl.cpp\n"
              << "  Functions:     write_rng_state(int) + read_rng_state(int)\n"
              << "  Scenario:      per-rank save/restore round trip (file guaranteed\n"
              << "                 to exist, since the failure path fcloses nullptr).\n"
              << "  Pass criteria: draws after the restore equal the draws that\n"
              << "                 immediately followed the save.\n";

    if (!randgsl_ensure_dir("DATA")) {
        std::cerr << "  !! Could not create DATA/ - skipping assertions\n";
        ADD_FAILURE() << "Could not create the DATA directory needed by write_rng_state(rank)";
        return;
    }

    const int rank = 5;
    const std::string src = "DATA/rng_state_" + std::to_string(rank);
    const std::string dst = "rng_state_" + std::to_string(rank);

    randgsl_init_rng(31337);
    for (int i = 0; i < 3; ++i)
        rand_gsl(); // advance so the state is non-trivial

    write_rng_state(rank);

    // read_rng_state(rank) looks in the current directory, not DATA/.
    ASSERT_TRUE(randgsl_copy_file(src.c_str(), dst.c_str()))
        << "Could not stage " << src << " as ./" << dst
        << "; refusing to call read_rng_state(rank) on a missing file";

    std::vector<double> expected;
    for (int i = 0; i < 4; ++i)
        expected.push_back(rand_gsl());

    for (int i = 0; i < 15; ++i)
        rand_gsl(); // move the generator far away from the saved point

    std::cerr << "  -> Calling read_rng_state(" << rank << ")\n";
    read_rng_state(rank);

    std::vector<double> restored;
    for (int i = 0; i < 4; ++i)
        restored.push_back(rand_gsl());

    for (size_t i = 0; i < expected.size(); ++i) {
        std::cerr << "     [" << i << "] expected " << expected[i]
                  << " vs restored " << restored[i] << "\n";
        EXPECT_DOUBLE_EQ(restored[i], expected[i])
            << "read_rng_state(rank) did not restore the saved state (draw " << i << ")";
    }

    // Leave the working directory tidy.
    remove(dst.c_str());
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers.  Each test_* helper is invoked from its own TEST so that
// a failure in one does not stop the others from running.
// -----------------------------------------------------------------------------
TEST(RandGsl, RandGslBasic) { test_randgsl_rand_gsl_basic(); }
TEST(RandGsl, RandGslUniformity) { test_randgsl_rand_gsl_uniformity(); }
TEST(RandGsl, RandGsl64) { test_randgsl_rand_gsl64(); }
TEST(RandGsl, SeedingSemantics) { test_randgsl_seeding_semantics(); }
TEST(RandGsl, GaussV) { test_randgsl_gaussv(); }
TEST(RandGsl, WriteRngState) { test_randgsl_write_rng_state(); }
TEST(RandGsl, WriteRngStateRank) { test_randgsl_write_rng_state_rank(); }
TEST(RandGsl, WriteRngStateSimItr) { test_randgsl_write_rng_state_simItr(); }
TEST(RandGsl, ReadRngStateRoundTrip) { test_randgsl_read_rng_state_roundtrip(); }
TEST(RandGsl, ReadRngStateMissingFile) { test_randgsl_read_rng_state_missing_file(); }
TEST(RandGsl, ReadRngStateRankRoundTrip) { test_randgsl_read_rng_state_rank_roundtrip(); }