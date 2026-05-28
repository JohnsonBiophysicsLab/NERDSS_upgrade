# Sanitizer And Static Analysis Guardrails

This page documents the Phase 2 local guardrails for serial NERDSS builds. The
commands intentionally do not change simulation behavior; they only change build
flags and analysis tooling.

## CMake Sanitizer Builds

Configure serial Debug builds with AddressSanitizer, UndefinedBehaviorSanitizer,
or both:

```sh
cmake -S . -B build/asan -DCMAKE_BUILD_TYPE=Debug -DNERDSS_ENABLE_ASAN=ON
cmake --build build/asan --target nerdss

cmake -S . -B build/ubsan -DCMAKE_BUILD_TYPE=Debug -DNERDSS_ENABLE_UBSAN=ON
cmake --build build/ubsan --target nerdss

cmake -S . -B build/sanitizers -DCMAKE_BUILD_TYPE=Debug \
  -DNERDSS_ENABLE_ASAN=ON \
  -DNERDSS_ENABLE_UBSAN=ON
cmake --build build/sanitizers --target nerdss
```

The helper below configures and builds all three variants under
`build/sanitizers/` by default:

```sh
tools/run_sanitizer_builds.sh
```

Useful environment variables:

- `NERDSS_SANITIZER_BUILD_ROOT`: override the helper build root.
- `NERDSS_BUILD_JOBS`: pass an explicit parallel job count to `cmake --build`.

## Make Sanitizer Builds

The Makefile keeps ordinary debug builds separate from sanitizer builds. Debug
objects now live under `obj/debug`, release objects under `obj/release`, and
sanitizer objects under `obj/asan`, `obj/ubsan`, or `obj/sanitizers`.

```sh
make asan
make ubsan
make sanitizers
```

The generated serial executables are:

- `bin/nerdss_asan`
- `bin/nerdss_ubsan`
- `bin/nerdss_sanitizers`

Sanitizer Make targets are serial-only. MPI sanitizer builds should be added
after the MPI build path has its own smoke coverage.

## Focused Static Analysis

Run `clang-tidy` on focused directories instead of the whole tree:

```sh
tools/run_static_analysis.sh src/math src/parser EXEs
```

If `.clang-tidy` is present, the script uses it. Until Agent E's style branch is
merged, the script falls back to a conservative high-signal check set:

```text
clang-analyzer-*, bugprone-*, performance-*, portability-*
```

Useful environment variables:

- `NERDSS_STATIC_ANALYSIS_BUILD_DIR`: override the compile database build dir.
- `NERDSS_CLANG_TIDY_CHECKS`: override the fallback check expression.

## Current Validation Notes

Validated locally on 2026-05-28 with AppleClang 16.0.0, CMake 4.3.3, and GSL
2.8 from Homebrew:

- `cmake -S . -B build/validate-asan -DCMAKE_BUILD_TYPE=Debug -DNERDSS_ENABLE_ASAN=ON`
  configures successfully.
- `cmake --build build/validate-asan --target nerdss --parallel 4` reaches
  `EXEs/nerdss.cpp` and then fails on existing master declarations listed in
  `docs/upgrade/static_analysis_backlog.md`.
- `cmake -S . -B build/validate-ubsan -DCMAKE_BUILD_TYPE=Debug -DNERDSS_ENABLE_UBSAN=ON`
  configures successfully.
- `make -n asan` shows the expected serial ASan compile and link flags and
  outputs `bin/nerdss_asan`.
- `tools/run_static_analysis.sh src/math` exits early because `clang-tidy` is
  not installed on the local PATH.
