# NERDSS Upgrade Implementation Plan

This plan splits the upgrade into assignable workstreams for dedicated agents.
The central rule is that core simulation behavior must remain unchanged unless a
separate, reviewed scientific change request explicitly approves it. Refactors,
format changes, diagnostics, and performance work must be validated against
baseline outputs before they land.

## Repository Observations

- Core engine is C++ with `src/` split by domain and public headers under
  `include/`.
- Current build entry points are `Makefile` and `CMakeLists.txt`; the Makefile
  supports `serial`, `mpi`, `debug`, and `profile`.
- Existing validation assets live mainly under `sample_inputs/VALIDATE_SUITE`
  and `run_code_tests`.
- Current outputs are mostly ad hoc `.dat`, `.xyz`, `.psf`, restart, and text
  files, with some JSON already present through `write_bonded_complex_json.cpp`
  and bundled `nlohmann_json`.
- Error handling currently mixes custom exception classes, direct `exit(...)`,
  and formatted `std::cerr` messages.
- Generated Doxygen HTML is committed under `docs/html`; future documentation
  should distinguish source docs from generated artifacts.

## Global Working Rules

1. Preserve behavior first. Every refactor branch starts by recording a baseline
   for at least one deterministic or fixed-seed example that exercises the code
   being touched.
2. Refactor in vertical slices. Prefer small, reviewable changes that keep the
   executable runnable over large rewrites.
3. Keep old input and output formats readable/writable until a documented
   deprecation window exists.
4. Separate mechanics from science. Formatting, naming, ownership, and data
   structure changes must not modify reaction probabilities, RNG streams,
   propagation, boundary behavior, or MPI partitioning unless explicitly scoped.
5. Each agent updates this document or a linked log under `docs/upgrade/` with
   decisions, validation runs, benchmarks, and unresolved risks.

## Phase 0: Baseline And Work Coordination

### Agent A: Repository Baseline And Branch Policy

Scope:
- Create a stable branch naming scheme, for example `codex/upgrade-*`.
- Record compiler, GSL, CMake, Make, OS, and optional MPI versions.
- Decide whether C++11 remains required or whether the project can move to
  C++14/C++17. The current CMake file uses C++11 while README badges mention
  C++ >=14, so this must be resolved before broad modernization.
- Add `docs/upgrade/decision_log.md` and `docs/upgrade/validation_log.md`.

Deliverables:
- Documented development environment matrix.
- Branch and PR template guidance.
- Decision on supported C++ standard.

Acceptance criteria:
- A new contributor can build the current serial executable from the docs.
- The team agrees on the minimum compiler and language standard.

Dependencies:
- None.

### Agent B: Baseline Build And Smoke Runner

Scope:
- Verify `make serial`, `make mpi` when MPI is available, and CMake serial build.
- Add a small script, for example `tools/run_smoke_tests.py`, that builds and
  runs a short sample input in an isolated temporary directory.
- Capture stdout, stderr, exit code, produced files, and runtime.

Deliverables:
- Repeatable smoke command for local and CI use.
- First baseline artifact bundle for selected samples.

Acceptance criteria:
- Smoke runner works on a clean checkout with documented dependencies.
- Build failures produce concise, actionable messages.

Dependencies:
- Agent A environment decision.

## Phase 1: Scientific Regression Harness

### Agent C: Deterministic And Fixed-Seed Validation

Scope:
- Select a minimal validation set from `sample_inputs/VALIDATE_SUITE`, starting
  with small systems such as `homoTrimer`, `hetTrimer`, `create_destroy`,
  `michaelis_menten`, `implicit_lipid`, and `sphere`.
- Standardize fixed seed invocation through the existing command parser.
- Define file-level comparisons for deterministic outputs and tolerant numeric
  comparisons for floating-point text outputs.
- Include restart validation, because `sample_inputs/VALIDATE_SUITE/README`
  explicitly requires checking both fresh starts and restart paths.

Deliverables:
- `tests/regression/` harness.
- Baseline manifests with file names, comparison mode, tolerances, and ignored
  volatile lines.
- Clear pass/fail report.

Acceptance criteria:
- The harness can compare current `master` against itself and pass.
- At least one fresh-run and one restart-run case are covered.

Dependencies:
- Agent B smoke runner.

### Agent D: Stochastic Coherence Validation

Scope:
- Define statistical checks for stochastic examples: means, ranges,
  distribution summaries, event counts, and monotonic invariants.
- Identify which outputs are scientifically meaningful for stochastic
  comparison and which are implementation artifacts.
- Run small ensembles at fixed seed sets to establish expected envelopes.

Deliverables:
- Stochastic validation spec.
- Ensemble runner and summary report.
- Initial thresholds reviewed by a domain expert.

Acceptance criteria:
- Stochastic validation detects large behavioral regressions without failing on
  normal random variation.

Dependencies:
- Agent C baseline format and runner conventions.

## Phase 2: Tooling, Style, And Static Safety

### Agent E: Google Style Tooling

Scope:
- Add `.clang-format` based on Google style, with local exceptions only where
  needed for readability or generated code.
- Add `.clang-tidy` checks for modernize, readability, bugprone, performance,
  and portability categories.
- Exclude generated docs, notebooks, binary assets, and vendored files.
- Create a formatting rollout plan by directory to avoid unreadable mega-diffs.

Deliverables:
- Style config files.
- `tools/format_changed_files.sh` or equivalent local command.
- First small formatting PR for one low-risk directory.

Acceptance criteria:
- Formatting can be checked in CI.
- Style changes are separated from behavior changes whenever possible.

Dependencies:
- Agent A C++ standard decision.
- Agent C regression harness for safety.

### Agent F: Static Analysis And Sanitizer Builds

Scope:
- Add AddressSanitizer/UndefinedBehaviorSanitizer debug targets through CMake
  and Make, preserving the existing debug target intent.
- Run `clang-tidy` on focused directories first.
- Inventory warnings by category and priority.

Deliverables:
- Sanitizer build documentation.
- Static analysis backlog.
- Initial fixes for high-confidence undefined behavior or memory safety issues.

Acceptance criteria:
- Sanitizer build compiles for serial mode.
- Known sanitizer failures are documented with reproduction commands.

Dependencies:
- Agent B build runner.
- Agent E style/tooling.

## Phase 3: Architecture Refactor Plan

### Agent G: Domain Model And SOLID Boundary Map

Scope:
- Map high-level responsibilities across parser, IO, system setup, reactions,
  trajectory, boundary conditions, math, MPI, and executable orchestration.
- Identify global state, hidden ownership, side effects, and large parameter
  lists that block testing.
- Propose target interfaces without changing algorithms yet.

Deliverables:
- `docs/upgrade/architecture_map.md`.
- Dependency diagram.
- Refactor priority list with risk levels.

Acceptance criteria:
- Each proposed extraction has a test or regression guard.
- The map identifies where scientific algorithms live versus orchestration and
  IO code.

Dependencies:
- Agent C regression harness.

### Agent H: Main Loop And Simulation Context Extraction

Scope:
- Extract setup, file initialization, simulation loop, output scheduling, and
  teardown from `EXEs/nerdss.cpp` into testable orchestration units.
- Introduce a `SimulationContext` or similarly named object only if it reduces
  global state and parameter threading.
- Keep RNG sequence and call order unchanged.

Deliverables:
- Smaller executable entry point.
- Testable setup and run functions.
- Documentation of preserved RNG behavior.

Acceptance criteria:
- Regression harness passes with fixed seeds.
- Output file names and default directories remain compatible.

Dependencies:
- Agent G architecture map.
- Agent C regression harness.

### Agent I: Parser Boundary Refactor

Scope:
- Separate lexical parsing, semantic validation, and model construction for
  `.inp`, `.mol`, restart, add-file, and BNGL-related parsing.
- Preserve old syntax and comments.
- Convert parser failures to structured errors with filename, line, key, and
  actionable message.

Deliverables:
- Parser module boundaries.
- Parser unit tests using small fixture files.
- Compatibility table for old input syntax.

Acceptance criteria:
- All existing sample inputs still parse.
- Invalid fixture files produce stable, clean messages.

Dependencies:
- Agent C validation harness.
- Agent K error model.

### Agent J: Reaction, Trajectory, And Boundary Hotspot Refactors

Scope:
- Refactor reaction and trajectory code only behind regression and benchmark
  guards.
- Prioritize repeated logic visible across box, sphere, compartment, implicit
  lipid, and MPI variants.
- Improve const-correctness, ownership, reserve patterns, and data locality
  without changing equations or RNG order.
- Split mechanical cleanup from performance changes.

Deliverables:
- Small PR series by subsystem:
  - Association and dissociation helpers.
  - Unimolecular and bimolecular reaction checks.
  - Boundary reflection helpers.
  - Trajectory propagation helpers.
  - Matrix/probability table lifetime management.
- Before/after benchmark notes for any performance-motivated change.

Acceptance criteria:
- Fixed-seed regression tests pass.
- Benchmarks show no unexplained slowdown.
- Domain expert review signs off on algorithm-sensitive changes.

Dependencies:
- Agent C regression harness.
- Agent M profiler baseline.
- Agent G architecture map.

## Phase 4: Standard Data IO With Backward Compatibility

### Agent K: Input Format Modernization

Scope:
- Design JSON schema for parameters, molecule templates, reactions,
  observables, restart metadata, and optional coordinate inputs.
- Provide converters from legacy `.inp` and `.mol` to JSON.
- Add parser dispatch that accepts legacy formats by default and JSON through
  explicit flags or file extension.
- Keep old inputs as first-class supported inputs during the upgrade.

Deliverables:
- `schemas/nerdss-input.schema.json`.
- Legacy-to-JSON converter.
- JSON parsing tests.
- Documentation with side-by-side legacy and JSON examples.

Acceptance criteria:
- A selected validation example can run from legacy files and equivalent JSON.
- Parsed internal model is equivalent between legacy and JSON paths.

Dependencies:
- Agent I parser boundary refactor.

### Agent L: Output Format Modernization

Scope:
- Define standard CSV outputs for time series currently written as `.dat`.
- Define JSON metadata manifest for each run: version, command, seed, input
  files, output files, schema versions, wall time, warnings, and git commit.
- Keep existing `.dat`, `.xyz`, `.psf`, PDB, and restart outputs by default
  unless a documented flag disables them.
- Add optional structured output directory layout under `DATA/`.

Deliverables:
- `schemas/nerdss-run-manifest.schema.json`.
- CSV writer helpers with stable headers.
- Run manifest writer.
- Compatibility documentation.

Acceptance criteria:
- Existing output files are still produced.
- CSV/JSON outputs can be consumed by Python without custom text parsing.
- Regression harness can compare both legacy and structured outputs.

Dependencies:
- Agent C regression harness.
- Agent H output scheduling boundaries.

## Phase 5: Error Handling And Diagnostics

### Agent M: Structured Error Model

Scope:
- Replace direct `exit(...)` calls incrementally with typed exceptions or
  status-return paths at module boundaries.
- Add error categories: input parse, file IO, invalid model, numerical failure,
  geometry failure, MPI failure, internal invariant violation.
- Ensure top-level executables convert errors to clean messages and stable exit
  codes.
- Include optional debug context without overwhelming normal users.

Deliverables:
- `include/error/` error type definitions.
- Top-level exception boundary in serial and MPI executables.
- Error message style guide.
- Tests for representative error cases.

Acceptance criteria:
- Invalid input fails without stack noise or segmentation fault.
- Internal invariant failures include enough context for developers.
- MPI rank context is preserved when relevant.

Dependencies:
- Agent I parser work for parse messages.
- Agent F sanitizer inventory.

### Agent N: Segfault Triage And Hardening

Scope:
- Run sanitizer builds and known problematic examples.
- Categorize crashes as input validation, memory ownership, vector/index bounds,
  null/dangling pointer, MPI synchronization, or numerical invalid value.
- Add assertions or checked access where the cost is acceptable, especially at
  module boundaries.
- Convert user-caused crashes into clean errors.

Deliverables:
- Crash reproduction catalog.
- Fixes for high-priority crashes.
- Regression tests for each fixed crash.

Acceptance criteria:
- Previously reproduced crashes either pass or fail cleanly.
- No new sanitizer regressions in covered smoke tests.

Dependencies:
- Agent F sanitizer builds.
- Agent M structured errors.

## Phase 6: Benchmarking And Profiling

### Agent O: Benchmark Harness

Scope:
- Select benchmark cases across small, medium, and representative large models.
- Measure wall time, CPU time, peak RSS where available, output size, and
  optional MPI scaling.
- Keep benchmark inputs separate from correctness validation when runtime is
  too high for CI.

Deliverables:
- `benchmarks/` runner.
- Baseline benchmark report in `docs/upgrade/benchmark_baseline.md`.
- Machine-readable results, preferably JSON or CSV.

Acceptance criteria:
- Benchmarks can be run locally with one command.
- Results include system metadata and git commit.

Dependencies:
- Agent B build runner.
- Agent C validation conventions.

### Agent P: Profiler And Hotspot Analysis

Scope:
- Use existing Makefile profiling intent and add documented alternatives such
  as gprof, perf, Instruments, or gperftools depending on platform.
- Profile at least one representative serial case and, if practical, one MPI
  case.
- Attribute time to reaction checks, propagation, overlap checks, IO, parsing,
  and MPI communication.

Deliverables:
- `docs/upgrade/profile_report.md`.
- Ranked hotspot list.
- Proposed optimization experiments with expected risk and validation needs.

Acceptance criteria:
- Heavy components are identified from measurements, not guesswork.
- Optimization tasks are linked to benchmark cases and regression guards.

Dependencies:
- Agent O benchmark harness.

## Phase 7: Automated Testing And GitHub Actions

### Agent Q: Unit Test Framework

Scope:
- Choose C++ test framework, for example GoogleTest or Catch2, considering
  dependency simplicity in GitHub Actions.
- Add unit tests for math primitives, parser helpers, format converters, error
  messages, and small pure functions first.
- Avoid testing scientific behavior only through mocks; keep numerical behavior
  anchored to regression tests.

Deliverables:
- `tests/unit/` structure.
- CMake test target with `ctest`.
- Initial unit tests for low-risk modules.

Acceptance criteria:
- Unit tests run locally and in CI.
- Tests do not require large sample assets unless marked integration.

Dependencies:
- Agent A C++ standard decision.
- Agent E tooling.

### Agent R: Integration And Regression CI

Scope:
- Add GitHub Actions workflows for formatting, build, unit tests, smoke tests,
  sanitizer build, and selected regression tests.
- Cache dependencies where appropriate.
- Include Linux first; add macOS once build times and dependencies are stable.
- Keep slow stochastic and benchmark suites manual or scheduled.

Deliverables:
- `.github/workflows/ci.yml`.
- Optional `.github/workflows/benchmark.yml`.
- CI documentation in README or contributor docs.

Acceptance criteria:
- PRs cannot merge when format, build, unit, or smoke tests fail.
- Regression failures report which files differ and how.

Dependencies:
- Agents B, C, E, F, Q.

## Phase 8: Documentation And Migration Support

### Agent S: User Documentation

Scope:
- Update install/build/run docs for serial and MPI.
- Document legacy and JSON/CSV input-output modes.
- Provide migration guide for old `.inp`, `.mol`, and `.dat` workflows.
- Keep examples short and runnable.

Deliverables:
- Updated README.
- `docs/input.md` updates.
- `docs/output.md` or equivalent.
- Migration examples.

Acceptance criteria:
- New users can run one serial example and inspect structured outputs.
- Existing users can keep legacy workflows without surprise breakage.

Dependencies:
- Agents K and L.

### Agent T: Developer Documentation And Backtracking

Scope:
- Maintain a chronological upgrade log.
- Record every meaningful refactor with motivation, files touched, validation
  command, result, and benchmark impact when applicable.
- Document known limitations and deferred work.

Deliverables:
- `docs/upgrade/change_log.md`.
- `docs/upgrade/known_risks.md`.
- Updated Doxygen/source documentation workflow.

Acceptance criteria:
- A reviewer can understand why each major change happened.
- Validation and benchmark evidence is linked from each major PR.

Dependencies:
- All agents.

## Suggested Agent Assignment Matrix

| Agent | Workstream | Primary paths | Output |
| --- | --- | --- | --- |
| A | Environment and policy | README, docs, build files | Environment and branch policy |
| B | Build and smoke runner | Makefile, CMakeLists, tools | Repeatable build/smoke command |
| C | Deterministic regression | sample_inputs, tests | Fixed-seed validation harness |
| D | Stochastic validation | sample_inputs, tests | Statistical coherence checks |
| E | Style tooling | repo root, src, include | clang-format/tidy configs |
| F | Static/sanitizer safety | CMake, Makefile, src | Sanitizer and warning backlog |
| G | Architecture map | docs, src, include | SOLID boundary plan |
| H | Main loop extraction | EXEs, src orchestration | Testable simulation runner |
| I | Parser refactor | src/parser, include/parser | Parser module and tests |
| J | Algorithm-adjacent refactor | src/reactions, trajectory, boundary | Guarded refactor PRs |
| K | Modern input | parser, schemas, tools | JSON schema and converter |
| L | Modern output | src/io, schemas | CSV/JSON outputs plus legacy |
| M | Error model | include/error, src/error, EXEs | Clean errors and exit codes |
| N | Crash hardening | src, tests | Segfault fixes and repro tests |
| O | Benchmarks | benchmarks, docs/upgrade | Benchmark harness |
| P | Profiling | benchmarks, docs/upgrade | Hotspot report |
| Q | Unit tests | tests/unit, CMake | Unit test framework |
| R | CI | .github/workflows | Automated checks |
| S | User docs | README, docs | Migration and usage docs |
| T | Backtracking docs | docs/upgrade | Change logs and risk register |

## Recommended Execution Order

1. Agents A, B, and C start first. No broad refactor should proceed until build,
   smoke, and deterministic validation exist.
2. Agents E and F add style/static/sanitizer guardrails once the build path is
   stable.
3. Agents G and O/P map architecture and performance hotspots in parallel.
4. Agents H, I, M, and N handle orchestration, parser, and diagnostics because
   they reduce risk for later work.
5. Agents K and L modernize input and output formats behind compatibility flags.
6. Agent J performs computationally adjacent refactors in small slices, guided
   by the profiler and guarded by regression tests.
7. Agents Q and R make the checks standard in local development and GitHub
   Actions.
8. Agents S and T continuously update user/developer documentation.

## PR Gate Checklist

Every PR in this upgrade should answer:

- What behavior is intentionally changed, if any?
- Which validation samples were run?
- Were deterministic outputs identical or within documented tolerance?
- Were stochastic outputs coherent under the documented envelope?
- Did runtime or memory change on the relevant benchmark?
- Were legacy inputs and outputs preserved?
- Are new errors cleaner than the old failure mode?
- Which docs or upgrade logs were updated?

## Initial Milestones

### Milestone 1: Guardrails

Target outcome:
- Serial build and smoke test are automated.
- Deterministic validation harness exists for at least two samples.
- Formatting and sanitizer configs exist.
- CI builds and runs smoke tests.

### Milestone 2: Compatibility Layer

Target outcome:
- Legacy parser behavior is covered by tests.
- JSON input schema and converter support one complete validation example.
- CSV/JSON structured outputs exist alongside legacy files.
- Clean error messages exist for common invalid inputs.

### Milestone 3: Safe Refactor

Target outcome:
- Main loop is split into testable units.
- Parser and IO boundaries are clearer.
- High-priority sanitizer/crash issues are resolved.
- Fixed-seed validation remains stable.

### Milestone 4: Performance Evidence

Target outcome:
- Benchmark harness and profile report identify top hotspots.
- At least one high-impact optimization experiment is completed or rejected with
  evidence.
- No unexplained quality or runtime regression is present.

### Milestone 5: Release-Ready Workflow

Target outcome:
- GitHub Actions cover format, build, unit, smoke, and selected regression tests.
- User migration docs are complete.
- Upgrade logs provide backtracking for major decisions and validation evidence.
