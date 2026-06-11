# C++ Unit Tests

Agent Q owns the Phase 7 unit test framework. The baseline CTest target uses a
small standard-library test executable and does not require Google Test.

## CTest Command

```bash
cmake -S . -B build-unit-tests
cmake --build build-unit-tests --target nerdss_unit_tests
ctest --test-dir build-unit-tests --output-on-failure
```

The first test target covers low-risk `Coord` and `Vector` helper behavior:
rounding, colinearity, magnitude, normalization, dot products, cross products,
projection, and angle calculation.

## Google Test Command

Google Test is only required for the explicit Makefile test target:

```bash
make unittest
```

If Google Test is not available through `pkg-config`, `make unittest` reports
the missing dependency and normal runtime builds remain available. On systems
without `pkg-config`, pass `GTEST_CFLAGS` and `GTEST_LIBS` explicitly.

Use `TEST_VERBOSE=1` for expanded passed-test details emitted by NERDSS tests:

```bash
make unittest TEST_VERBOSE=1
```

Additional Google Test flags can be forwarded with `GTEST_FLAGS`, for example:

```bash
make unittest GTEST_FLAGS="--gtest_filter=NormFunctionTest.* --gtest_print_time=1"
```

Generated build directories such as `build-unit-tests/` should not be
committed.
