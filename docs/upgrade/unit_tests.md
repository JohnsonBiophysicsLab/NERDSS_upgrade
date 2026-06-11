# C++ Unit Tests

Agent Q owns the Phase 7 unit test framework. For the current Google Test
cleanup, the supported path is Makefile-only.

## Scope Warning

Do not add or extend CMake/CTest integration for Google Test in this slice.
That work is out of scope for now. Use the Makefile test target below, and only
reopen CMake work after a separate reviewed request explicitly asks for it.

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
