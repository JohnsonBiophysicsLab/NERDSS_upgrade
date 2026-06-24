// Main entry point for the norm_function GoogleTest suite.
// Defines the globals required by linked translation units (rand_gsl.o,
// find_which_reaction.o, etc.) that declare them extern, then hands control
// to the GTest runner.

#include <gsl/gsl_rng.h>
#include <gtest/gtest.h>

gsl_rng* r = nullptr;
unsigned long totMatches = 0;

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
