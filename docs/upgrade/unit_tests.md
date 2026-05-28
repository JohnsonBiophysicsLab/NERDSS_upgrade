# C++ Unit Tests

Agent Q owns the Phase 7 unit test framework. The initial framework uses CTest
with a small standard-library test executable, avoiding network-fetched test
dependencies.

## Command

```bash
cmake -S . -B build-unit-tests
cmake --build build-unit-tests --target nerdss_unit_tests
ctest --test-dir build-unit-tests --output-on-failure
```

The first test target covers low-risk `Coord` and `Vector` helper behavior:
rounding, colinearity, magnitude, normalization, dot products, cross products,
projection, and angle calculation.

Generated build directories such as `build-unit-tests/` should not be
committed.
