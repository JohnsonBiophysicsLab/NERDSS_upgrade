// Global definitions that normally live in the EXE main (EXEs/nerdss.cpp).
// The project object files reference these as `extern`, so any test that links
// against the full object set needs them defined exactly once here.

#include <gsl/gsl_rng.h>

gsl_rng* r;             // global GSL random generator (see math/rand_gsl.hpp)
long long randNum = 0;  // random-number call counter
unsigned long totMatches = 0; // reaction-match counter (see classes/class_Rxns.hpp)
