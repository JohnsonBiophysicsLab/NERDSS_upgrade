# Phase 2 Style Tooling

This branch adds the initial Google-style formatting and clang-tidy tooling for
NERDSS. It intentionally does not bulk-format `src/` or `include/`; behavior
changes and mechanical style changes should stay in separate pull requests.

## Formatting

- Root `.clang-format` is based on Google C++ style with C++11 parsing and a
  100-column limit.
- Use `tools/format_changed_files.sh` to format only C/C++ files changed by the
  current branch:

  ```sh
  tools/format_changed_files.sh
  tools/format_changed_files.sh --check
  ```

- For a staged rollout by directory, use small PRs in this order:
  1. `EXEs/` and focused helper-only files.
  2. Low-risk leaf directories such as `src/debug/`, `src/error/`, and matching
     headers.
  3. Domain directories only after regression coverage is available for the
     touched behavior.

## Static analysis

- Root `.clang-tidy` enables practical analyzer, `bugprone`, `performance`,
  `portability`, selected `modernize`, and selected `readability` checks.
- The config assumes C++11 until the project-wide standard decision changes.
- Run clang-tidy on focused compile-command entries rather than the full tree
  while the upgrade is in progress.

## Exclusions

Do not format or run clang-tidy over generated, external, or non-source assets:

- `docs/html/**` generated Doxygen output.
- `third_party/**` vendored dependencies.
- `**/*.ipynb` notebooks.
- Binary/reference assets such as `*.pdf`, `*.png`, `*.jpg`, `*.svg`, `*.eps`,
  and `*.docx`.
- Build and IDE output directories such as `build/`, `cmake-build-*`, `bin/`,
  `obj/`, `.idea/`, and `.vscode/`.

The `.clang-format-ignore` file records these formatting exclusions for
clang-format versions that support it; the helper script applies the same
exclusion policy explicitly.
