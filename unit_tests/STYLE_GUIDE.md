# NERDSS Unit Test Style Guide

How to write a unit test that matches the 192 test files already in
[unit_tests/src/](src/). Every rule below is derived from those files; where two
patterns exist in the tree, the one documented here is the one to follow going
forward.

The short version: **one test file per source file, a Doxygen header that
explains the function under test, helpers in an anonymous namespace, scenarios
in prefixed `void test_*()` functions that narrate themselves to `stderr`,
`EXPECT_*` assertions with explanatory messages, and a block of one-line
`TEST()` wrappers at the bottom.**

---

## 1. Where things live and how they build

| Thing | Path |
| --- | --- |
| Test sources | `unit_tests/src/test_<source_file_name>.cpp` |
| Entry point (do not duplicate) | `unit_tests/src/gtest_main.cpp` |
| Objects | `unit_tests/obj/` |
| Binary | `unit_tests/bin/gtest_main` |

```sh
make unittest                  # from the repository root
./unit_tests/bin/gtest_main    # run everything

# Run one file's tests, or one scenario:
./unit_tests/bin/gtest_main --gtest_filter='AssociateImplicitLipidBox.*'
./unit_tests/bin/gtest_main --gtest_filter='ClassCoordTest.Roundv'
./unit_tests/bin/gtest_main --gtest_list_tests
```

`UTSRCS = $(wildcard $(UTSDIR)/*.cpp)` in the [Makefile](../Makefile), so a new
file in `unit_tests/src/` is picked up automatically — there is no list to edit.

### Constraints that follow from the build

These are not stylistic. Breaking one of them breaks the whole suite.

- **C++17** (`CXXFLAGS = -std=c++17`). No newer language or library features.
  Note this is looser than the C++11 target for `src/` in
  [docs/upgrade/contributor_guidance.md](../docs/upgrade/contributor_guidance.md);
  it applies to test code only.
- **One binary, one namespace for every test file.** Every free function name
  you define must be unique across all 192 files. Prefix everything
  (see §4).
- **Never define `main`.** `gtest_main.cpp` provides it.
- **Never define `r` or `totMatches`.** `gtest_main.cpp` defines them; other
  translation units declare them `extern`. A second definition is a link error
  for the entire suite. Use `extern gsl_rng* r;` or include
  `"math/rand_gsl.hpp"`.
- **`include/` is on the include path.** Use `#include "classes/class_Vector.hpp"`.
  Never `../include/...`, never a bare `#include "class_Vector.hpp"`.
- **gtest comes from pkg-config** → `#include <gtest/gtest.h>` is correct.
- **GSL is linked in**, so `<gsl/gsl_rng.h>` and friends are available.

---

## 2. The canonical file skeleton

```cpp
/*! \file test_<source_name>.cpp
 *
 * ### Unit test for src/<dir>/<source_name>.cpp
 *
 * The function under test is
 *
 *     <exact signature, wrapped and indented>
 *
 * <What the routine does, in prose. Name the distinct code paths and say
 *  which ones the tests below reach.>
 *
 * <Any precondition the tests have to satisfy, e.g. "Every test builds a
 *  complete two-molecule system so that no member read by the routine is
 *  left uninitialized.">
 */

#include "<dir>/<header_declaring_the_function>.hpp"   // project headers first

#include <gtest/gtest.h>                               // then gtest

#include <cmath>                                       // then std, alphabetical
#include <iostream>
#include <vector>

namespace {

// -----------------------------------------------------------------------------
// Fixed indices / sizes / tolerances shared by all tests in this file.
// -----------------------------------------------------------------------------
constexpr double kPfxTol = 1e-9;
constexpr int kPfxProtMolIdx = 0;   //!< moleculeList index of the protein

/*! \brief <One-line description of the helper.>
 *
 * \param[in]  <name> <meaning>
 * \param[out] <name> <meaning>
 */
ReturnType pfx_helper(Args...)
{
    ...
}

} // namespace

// -----------------------------------------------------------------------------
// Test 1: <one-line description of the scenario>.
// -----------------------------------------------------------------------------
void test_pfx_first_scenario()
{
    std::cerr << "\n[TEST] test_pfx_first_scenario\n"
              << "  Source file:   src/<dir>/<source_name>.cpp\n"
              << "  Function:      <function_name>()\n"
              << "  Scenario:      <what is being set up>\n"
              << "  Pass criteria: <what must hold afterwards>\n";

    // Arrange: build a fully initialized system.
    ...

    std::cerr << "  Calling <function_name>()...\n";
    function_under_test(...);

    // Assert, with a message on every check.
    EXPECT_NEAR(result.x, 1.0, kPfxTol) << "x must be unchanged by the no-op branch";
}

// -----------------------------------------------------------------------------
// GoogleTest wrappers. Each scenario is reported separately, and because the
// helpers use EXPECT_* (not ASSERT_*) every test runs to completion even when
// one check fails.
// -----------------------------------------------------------------------------
TEST(PfxSuiteName, FirstScenario) { test_pfx_first_scenario(); }
TEST(PfxSuiteName, SecondScenario) { test_pfx_second_scenario(); }
```

Reference implementations to read before writing a new file:
[test_class_Coord.cpp](src/test_class_Coord.cpp) (simple, class-based),
[test_reflect_traj_check_span.cpp](src/test_reflect_traj_check_span.cpp)
(free function with helpers), and
[test_associate_ImplicitLipid_box.cpp](src/test_associate_ImplicitLipid_box.cpp)
(large system setup).

---

## 3. The file header comment

Always a Doxygen `/*! \file ... */` block, always the first thing in the file.
178 of the 192 existing files have one. It must contain, in order:

1. `\file test_<name>.cpp`
2. `### Unit test for src/<dir>/<name>.cpp`
3. The **exact signature** of the function under test, indented as a code block.
4. **What the function does**, written from reading the implementation — not
   guessed from the name. Enumerate the branches:
   `1. Both complexes already on the membrane -> ... 2. The binder is in
   solution -> ...`.
5. **What the tests verify**, and anything deliberately *not* covered (an
   `exit()` path, an unreachable branch) with the reason.
6. Any environmental requirement — e.g. `test_write_crds.cpp` documents that it
   must create the relative `out/` directory first or the `ofstream` silently
   fails.

This block is the one place a reader learns why the assertions are what they
are. Write it last, after you know what the code actually does, but put it
first in the file.

---

## 4. Naming

Everything in a test file is linked into one binary with everything else, so
names collide unless they are prefixed. Pick a short lowercase abbreviation of
the source file name and use it everywhere in that file.

| Element | Rule | Example |
| --- | --- | --- |
| File | `test_` + exact source file name | `test_break_interaction.cpp` |
| Prefix | short abbreviation of the source name | `bi_`, `ailbox_`, `rtcs_`, `wc_` |
| Scenario function | `void test_<prefix>_<snake_case_scenario>()` | `test_bi_dimer_full_dissociation()` |
| Helper | `<prefix>_<verb>_<noun>()` | `ailbox_build_system()`, `rtcs_make_molecule()` |
| Constant | `k` + CamelCase with prefix | `kAilBindRadius`, `kRotateFnTol` |
| TEST suite | CamelCase of the source name | `BreakInteractionTest`, `AssociateImplicitLipidBox` |
| TEST case | CamelCase of the scenario | `DimerFullDissociation` |
| Death-test suite | suite name **must** end in `DeathTest` | `ErrorDeathTest` |

The `DeathTest` suffix is a GoogleTest requirement, not a preference: it makes
gtest run those suites first so forked children do not inherit a dirty heap.

Suite names are reasonably consistent but not uniform (`ClassCoordTest` vs
`AssociateSphere`). Either is acceptable; be consistent within a file, and use
exactly one suite name per file so `--gtest_filter='Suite.*'` selects that file.

---

## 5. Structure within the file

### Helpers go in an anonymous namespace; scenarios do not

174 of the 192 files wrap helpers, constants, and fixture structs in
`namespace { ... } // namespace`. The `void test_*()` scenario functions sit
*outside* it, above the `TEST()` wrappers. Keep that split: the namespace is
for things that must not leak, and the scenario functions read better at file
scope next to the wrappers that call them.

### Separator banners

Every logical section is introduced by a 79-column rule comment:

```cpp
// -----------------------------------------------------------------------------
// Test 2: the payload line is formatted as "Error: <message>!!!".
// -----------------------------------------------------------------------------
```

Use them for: the helper block, each scenario function, and the wrapper block.

### Do not use `TEST_F` fixtures

Zero existing files use `TEST_F` or derive from `::testing::Test`. When several
scenarios need the same world, build a plain struct plus a builder function
instead — see `AilSystem` / `ailbox_build_system()` in
[test_associate_ImplicitLipid_box.cpp](src/test_associate_ImplicitLipid_box.cpp):

```cpp
struct AilSystem {
    Parameters params {};
    Membrane membrane {};
    std::vector<MolTemplate> molTemplateList {};
    ...
};

void ailbox_build_system(AilSystem& sys, bool proteinComplexOnSurface,
    bool proteinIsPoint, bool ilIsFirstReactant, const Coord& proteinCom);
```

Give the builder `bool`/enum knobs for the branches you want to select, and
call it with named-argument comments at the call site so the booleans stay
readable:

```cpp
ailbox_build_system(sys, /*onSurface*/ true, /*isPoint*/ false,
                    /*ilFirst*/ false, protStart);
```

### Wrappers at the bottom, one line each

```cpp
TEST(AssociateImplicitLipidBox, TwoDBindingMol2IsImplicitLipid) { test_ailbox_2d_mol2_is_implicit_lipid(); }
```

No logic in the wrapper. This keeps the scenario functions ordinary C++ that
you can read, reorder, and call from each other if ever needed.

---

## 6. Console output

These tests are meant to be readable as a running log, not just a pass/fail
tally. Everything goes to `std::cerr` (unbuffered, interleaves correctly with
gtest's own output).

**Every scenario opens with a four-field banner:**

```cpp
std::cerr << "\n[TEST] test_bi_dimer_full_dissociation\n"
          << "  Source file:   src/reactions/break_interaction.cpp\n"
          << "  Function:      break_interaction()\n"
          << "  Scenario:      bound dimer A(a!1).B(b!1) dissociates; the two\n"
          << "                 molecules are not connected by any other bond.\n"
          << "  Pass criteria: returns false (not cancelled), breakLinkComplex\n"
          << "                 false, both interfaces freed.\n";
```

Keep the labels aligned exactly as above and continue wrapped lines under the
text column. `Pass criteria` is not optional — it is what lets a reader judge
whether a failure is a real regression or a bad test.

**Inside the scenario:**

```cpp
std::cerr << "  -> Testing roundv() with negative value -1.234567\n";  // step
std::cerr << "  Calling reflect_traj_check_span...\n";                 // the call
std::cerr << "  Final radial distance after reflection = " << radial
          << " (sphere radius " << membraneObject.sphereR << ")\n";    // computed value
std::cerr << "  Checked: exit status == 1 and message reached stderr.\n";
```

Print any value you compute to decide an assertion. When a test fails six
months from now, that number is the first thing someone needs.

Simple per-function test files (`test_class_Coord.cpp`) may use the lighter
form — `[<source>.cpp] Testing <fn>()...` followed by `  -> ` step lines —
instead of the full banner. Use the full banner whenever the scenario involves
setup that a reader would otherwise have to reconstruct.

---

## 7. Assertions

### `EXPECT_*`, not `ASSERT_*`

A fatal assertion aborts the current function, which in this suite means the
remaining checks for that scenario are silently skipped. Default to `EXPECT_*`
— the existing files use it 7300 times against 466 `ASSERT_*`.

Use `ASSERT_*` only where continuing would be undefined behavior rather than
merely uninformative, and say so:

```cpp
// The routine must have produced a second complex; indexing it below would
// read out of bounds if it did not.
ASSERT_EQ(complexList.size(), 2u) << "dissociation must create a new complex";
```

### Picking the macro

| Comparison | Macro |
| --- | --- |
| Integers, enums, sizes, strings, exact counters | `EXPECT_EQ` |
| Doubles produced by an exactly-specified operation (an assignment, a documented rounding rule) | `EXPECT_DOUBLE_EQ` |
| Any double that went through arithmetic | `EXPECT_NEAR(actual, expected, kPfxTol)` |
| Booleans | `EXPECT_TRUE` / `EXPECT_FALSE` |
| Bounds, "stayed inside the box" | `EXPECT_LE` / `EXPECT_GE` |
| A call that terminates the process | `EXPECT_EXIT` (see §9) |

Never write a bare numeric literal as a tolerance. Declare one `constexpr
double k<Prefix>Tol` near the top and use it everywhere in the file —
`kAssocBoxTol = 1e-9`, `kMtxfnTol = 1e-12`, `kDnTol = 1e-6`. One named value
per file makes the precision claim visible and adjustable.

### Every assertion carries a message

```cpp
EXPECT_DOUBLE_EQ(c1.x, 1.0) << "x should be 1.0";
EXPECT_EQ(protein.interfaceList[0].interaction.conjBackRxn, -1)
    << "irreversible reactions must not store a back-reaction index";
```

The message states the *expected behavior in words*, not a restatement of the
operands. Omit it only for a check that is self-evident from the line above it
in a tight group.

### Assert against the implementation, not against intuition

Read the source before choosing an expected value. `Coord::operator==` compares
components *after* `roundv()`, which truncates to four decimals with a
sign-dependent rule; picking your own epsilon instead produces an assertion
that is wrong about the code. If a function's behavior is surprising, encode
the surprising behavior and explain it in a comment.

---

## 8. Random number generation

Never call `srand_gsl()`. Initialize the shared global generator with this
exact idiom, guarded so repeated calls from multiple scenarios are harmless:

```cpp
/*! \brief Initialise the global GSL generator once. */
void pfx_init_rng()
{
    if (r == nullptr) {
        const gsl_rng_type* T;
        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 42);
    }
}
```

`r` is defined in `gtest_main.cpp` and declared `extern` in
`math/rand_gsl.hpp`. Include that header, or write `extern gsl_rng* r;`
yourself — do not define it.

Call `pfx_init_rng()` at the top of every scenario that reaches code which may
draw a random number, even indirectly. Seed 42, always: tests must be
deterministic across runs and across `--gtest_filter` selections. Do not assert
on specific random draws; assert on invariants the draw cannot violate.

---

## 9. Not crashing the suite

All 192 test files share one process. A segfault or `exit()` anywhere takes
down every remaining test, so these rules matter more than usual.

**Fully initialize everything you construct.** NERDSS constructors and member
functions routinely index into vectors without bounds checks and read members
the caller was expected to fill in. An under-initialized `Molecule`, `Complex`,
`Parameters` or `MolTemplate` segfaults the binary. Read the implementation for
the preconditions and satisfy all of them — this is why the larger tests use a
single `*_build_system()` helper rather than ad-hoc setup per scenario.

**Never reach an `exit()` or `abort()` path inline.** Several constructors bail
out that way on bad input (`Coord(std::vector<double>)` calls `exit(1)` when
handed more than three values). If a path like that exists, either note it in a
comment instead of exercising it, or test it properly as a death test:

```cpp
// ExitedWithCode(1) fails the test if the child crashed by signal or returned
// any status other than 1. The regex argument additionally requires the given
// text to appear on the child's stderr.
EXPECT_EXIT(error(std::string { "fatal condition detected" }),
    ::testing::ExitedWithCode(1),
    "Error: fatal condition detected")
    << "error(std::string) must call exit(1) and report the message";
```

The third argument is a POSIX extended regex matched against the child's
stderr. If your expected text contains regex metacharacters, either escape them
or pick a different substring — and say in a comment that you checked
(`'!' and ':' are not metacharacters, so this is a plain substring search`).
Remember the suite name must end in `DeathTest`.

---

## 10. Tests that touch the filesystem

Several routines under test write files. Such a test must:

1. **Create the output directory first.** `write_crds()` opens `out/<name>.dat`
   relative to the working directory; if `out/` is missing the `ofstream`
   construction fails silently and nothing is written — the test would pass
   vacuously. Call `mkdir("out", 0755)` (ignore `EEXIST`) in a helper.
2. **Use a filename unique to the test file**, since all tests share one
   working directory.
3. **Remove what it created**, with `std::remove(name.c_str())` at the end of
   the scenario.
4. **Assert on contents, not just existence** — read the file back with an
   `ifstream` and check the specific fields the routine is documented to write.

---

## 11. Formatting

Test code follows the layout in
[docs/devel.md](../docs/devel.md) — WebKit-ish, which is what the existing
tests are written in:

- 4-space indent.
- Opening brace **on its own line for function definitions**, attached at the
  end of the line everywhere else (`if`, `for`, `struct`, lambdas).
- Column limit 120.
- `Type* ptr` and `Type& ref` (pointer/reference bound to the type).
- Braced initialization where it reads well: `Coord { 0.0, 0.0, 0.0 }`,
  `long tot { 0 }`.

Include order, in three blocks separated by blank lines:

1. Project headers (`"classes/...hpp"`, `"reactions/...hpp"`), alphabetical.
2. `<gtest/gtest.h>`.
3. Standard library headers, alphabetical.

Include only what you use, and add a trailing comment on any include whose
purpose is not obvious: `#include "math/rand_gsl.hpp"  // extern gsl_rng* r`.

> **Note:** the repository `.clang-format` specifies Google style with 2-space
> indent, which does *not* match the existing tests. Do not run clang-format
> over `unit_tests/src/` — it would reformat 192 files and bury real diffs.
> Match the surrounding code by hand.

---

## 12. Comment density

These files are commented far more heavily than typical production code, and
that is deliberate: a test is the documentation of a behavior. Expect a comment
roughly every 3–5 lines of test body. Specifically:

- A `/*! \brief ... */` Doxygen block on every helper, with `\param[in]` /
  `\param[out]` / `\return` where there is anything to say.
- A comment above each setup block explaining *why* the value is what it is:
  `// scaleMaxDisplace is set huge so that measure_complex_displacement() can
  never reject the move.`
- A comment on any value chosen to select a code path:
  `prot.checkOverlap = false; // keeps the system-wide overlap check trivial`
- `//!<` trailing comments on constants and struct members.
- A comment wherever the test accommodates a quirk of the implementation, so
  nobody later "fixes" the test.

---

## 13. Coverage expectations for a new file

- Test **every function** defined in the source file, including every overload
  and every constructor.
- Cover **each branch** you can reach safely, and name in the file header the
  ones you could not and why.
- Include degenerate inputs that the code is expected to survive: empty vectors,
  empty `memberList`, zero-length trajectories, identical input coordinates.
- Include at least one invariant check per scenario that is independent of the
  exact numbers — "the complex ended up inside the sphere", "the three event
  histograms sum to one event" — so the test survives harmless numerical
  changes.
- Where a routine has symmetric arguments, test the reversed order. Reactant
  ordering bugs are a recurring class of defect here (see
  `BreakInteractionTest.DimerReversedReactantOrder`).

---

## 14. Checklist before committing

- [ ] File is `unit_tests/src/test_<exact source name>.cpp`.
- [ ] `/*! \file */` header names the source file, gives the signature, and
      describes the branches.
- [ ] Every free function and constant carries the file's unique prefix.
- [ ] Helpers and constants are inside `namespace { ... } // namespace`.
- [ ] No `main`, no definition of `r` or `totMatches`.
- [ ] Project includes use paths relative to `include/`.
- [ ] RNG initialized with the `gsl_rng_default` / seed-42 idiom, not
      `srand_gsl()`.
- [ ] `EXPECT_*` throughout; each `ASSERT_*` justified in a comment.
- [ ] Floating-point comparisons use a named `constexpr double k*Tol`.
- [ ] Each scenario prints the `[TEST]` banner with `Pass criteria`.
- [ ] Any file written is created in an existing directory and removed after.
- [ ] `make unittest` builds clean, and
      `./unit_tests/bin/gtest_main --gtest_filter='YourSuite.*'` passes.
- [ ] The full suite still passes — a link error or a crash in your file breaks
      everyone else's tests.

---

## 15. Generating a first draft

[generate_multi_unit_test_gtest.py](generate_multi_unit_test_gtest.py) produces
a starting file from a source path, using the whole of `include/` and
`src/classes/` as context plus `test_class_Coord.cpp` and
`test_reflect_traj_check_span.cpp` as style references.

```sh
cd unit_tests
python generate_multi_unit_test_gtest.py ../todo_files.txt          # skips existing
python generate_multi_unit_test_gtest.py ../todo_files.txt --force  # regenerate
```

Output is a **draft**. It has historically needed hand-checking for
under-initialized objects that segfault, assertions asserted against the name
of a function rather than its behavior, and branches that cannot actually be
reached with the setup provided. Build and run it before committing, and treat
the rules above as the acceptance criteria.

If you change the conventions in this guide, update `SYSTEM_PROMPT` /
`BUILD_CONVENTIONS` in that script too, or the next generated batch will drift
back.
