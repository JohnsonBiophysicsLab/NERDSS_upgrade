# Writing a Unit Test for NERDSS

A short, practical manual. It assumes you know NERDSS — molecules, complexes,
reactions, the `Parameters` object — but not GoogleTest. Read this once, write
your first test, then use [STYLE_GUIDE.md](STYLE_GUIDE.md) as the reference for
details.

---

## 1. What you are actually doing

A NERDSS unit test calls **one function from `src/`** with inputs you construct
by hand, and checks that what comes back is what the implementation promises.
No input files, no simulation loop, no `nerdss` executable — just the function.

All the tests live in [unit_tests/src/](src/), one file per source file, and
they are compiled together into **a single program**, `unit_tests/bin/gtest_main`.
That one-program detail drives most of the rules in §6: your test shares a
process with 191 others, so a crash in your file takes all of them down.

The framework is **GoogleTest** (gtest). You need about four things from it,
covered in §5.

---

## 2. First, run the suite that already exists

Before writing anything, confirm you can build and run:

```sh
cd <nerdss root>
make unittest                     # builds unit_tests/bin/gtest_main
./unit_tests/bin/gtest_main       # runs all 1443 tests, ~3 minutes
```

The first build is slow because it compiles all of `src/` plus all the tests.
After that it is incremental — adding one test file recompiles one file and
relinks.

To run a subset (you will do this constantly):

```sh
./unit_tests/bin/gtest_main --gtest_filter='AreSameAngleTest.*'
./unit_tests/bin/gtest_main --gtest_filter='AreSameAngleTest.ToleranceBoundary'
./unit_tests/bin/gtest_main --gtest_list_tests        # see what exists
```

---

## 3. The procedure, in seven steps

1. **Pick a function.** Start small: one function, few arguments, no `Molecule`
   or `Complex` setup. [../todo_files.txt](../todo_files.txt) lists source files
   that still have no test.
2. **Read its implementation line by line.** This is the real work. You are
   writing down what the code *does*, not what its name suggests. Note every
   branch, every member of every struct it reads, and anything surprising.
3. **Create `unit_tests/src/test_<source file name>.cpp`.** Exact same name as
   the source file, with `test_` in front. Nothing to register — the Makefile
   globs the directory.
4. **Write the file header comment** stating the signature, what the function
   does, and which branches you will cover.
5. **Write one `void test_<prefix>_<scenario>()` function per scenario**, each
   printing a banner saying what it is testing, then setting up inputs, calling
   the function, and checking the result with `EXPECT_*`.
6. **Add one-line `TEST()` wrappers at the bottom** so gtest can find and run
   them.
7. **`make unittest`, run your filter, then run the whole suite** to confirm you
   did not break anyone else.

---

## 4. A complete worked example

Read [src/test_areSameAngle.cpp](src/test_areSameAngle.cpp) alongside this
section. It is a good first model: the function under test is two lines, so
nothing is hidden, and the test still shows every habit the larger files use.

The function,
[src/reactions/areSameAngle.cpp](../src/reactions/areSameAngle.cpp):

```cpp
bool areSameAngle(double ang1, double ang2)
{
//    return std::abs(ang1) == std::abs(ang2);
    return std::abs(ang1 - ang2) < 1E-4;
}
```

**Step 2, reading it, is where the test gets its value.** Four things fall out
that the name alone would not tell you — and each one becomes a test:

- The comparison is **strict `<`**, so a difference of *exactly* `1E-4` means
  "different". That boundary is the single most likely thing to change by
  accident.
- The tolerance is **absolute, not relative**, so it behaves the same near `0`
  as near `pi`.
- The commented-out line shows the function *used to* compare magnitudes. It no
  longer does, so `pi` and `-pi` are different angles — and there is no `2*pi`
  wrap-around either.
- `std::abs(NaN - NaN)` is `NaN`, and `NaN < 1E-4` is false, so **NaN is never
  "the same angle" as anything, including itself.** This matters because the
  association code's default angles are quiet NaN.

Those observations are written into the file header, which is the first thing
in the file:

```cpp
/*! \file test_areSameAngle.cpp
 *
 * ### Unit test for src/reactions/areSameAngle.cpp
 *
 *     bool areSameAngle(double ang1, double ang2)
 *
 * The implementation is a tolerance-based floating point comparison:
 *
 *     return std::abs(ang1 - ang2) < 1E-4;
 *
 * ...
 *
 * Note that the function does NOT compare magnitudes (the commented-out
 * legacy behaviour), so an angle and its negation are only "the same" when
 * the angle is (nearly) zero. There is also no 2*pi wrap-around handling,
 * so pi and -pi are reported as different. Both of those properties are
 * pinned down explicitly by the tests below so a future behaviour change is
 * caught.
 */
```

That last sentence is the point of a unit test in one line: **it pins down
current behaviour so that changing it has to be deliberate.**

### The parts of the file

**Includes**, in three blocks — project, then gtest, then standard library:

```cpp
#include "reactions/association/association.hpp"

#include <gtest/gtest.h>

#include <cmath>
#include <iostream>
```

**A helper in an unnamed namespace.** Seven scenarios all do "call the function
with two angles, compare against an expectation, say what happened", so that is
factored out once. The `asa_` prefix is what keeps the name from colliding with
the other 190-odd test files:

```cpp
namespace {

/*! \brief Call areSameAngle() and check the result against an expectation.
 *
 * \param[in] ang1     first angle passed to areSameAngle()
 * \param[in] ang2     second angle passed to areSameAngle()
 * \param[in] expected the value areSameAngle() is required to return
 * \param[in] why      human readable description of the pass criterion
 */
void asa_check(double ang1, double ang2, bool expected, const char* why)
{
    const bool actual = areSameAngle(ang1, ang2);
    std::cerr << "    areSameAngle(" << ang1 << ", " << ang2 << ") -> "
              << std::boolalpha << actual
              << "  | expected " << expected
              << "  | |diff| = " << std::abs(ang1 - ang2)
              << "  | " << why << '\n';

    EXPECT_EQ(actual, expected)
        << "areSameAngle(" << ang1 << ", " << ang2 << ") returned " << actual
        << " but " << expected << " was required (" << why << ")";
}

} // namespace
```

**One function per scenario**, each opening with the four-field banner:

```cpp
// -----------------------------------------------------------------------------
// Test 4: the boundary itself. The implementation uses a strict '<', so a
//         difference of exactly 1E-4 must be reported as different.
// -----------------------------------------------------------------------------
void asa_test_tolerance_boundary()
{
    std::cerr << "\n[TEST] asa_test_tolerance_boundary\n"
              << "  Source file:   src/reactions/areSameAngle.cpp\n"
              << "  Function:      areSameAngle(double, double)\n"
              << "  Scenario:      difference is exactly the tolerance 1E-4.\n"
              << "  Pass criteria: strict '<' comparison means false.\n";

    // 0.0 and 1e-4 are both exactly representable operations here, so the
    // subtraction yields exactly 1e-4 and the strict comparison fails.
    asa_check(0.0, 1e-4, false, "exactly 1e-4 apart -> not within strict tolerance");
    asa_check(1e-4, 0.0, false, "same boundary case with arguments swapped");
    asa_check(0.0, -1e-4, false, "exactly -1e-4 apart -> also outside");
}
```

Note the comment justifying *why* these particular numbers are safe: `0.0` and
`1e-4` subtract exactly, so the test is not at the mercy of rounding. Choosing
inputs that make the arithmetic exact is a habit worth copying.

**One scenario that checks a property rather than a value.** Scenario 5 sweeps
nine angles and asserts only that the function is symmetric and reflexive — no
hardcoded expectations at all, so it keeps working if the tolerance is ever
retuned:

```cpp
const double angles[] = { -M_PI, -1.5, -1e-5, 0.0, 1e-5, 0.5, M_PI_2, M_PI, 10.0 };

for (int i = 0; i < numAngles; ++i) {
    EXPECT_TRUE(areSameAngle(angles[i], angles[i]))
        << "angle " << angles[i] << " must equal itself";

    for (int j = 0; j < numAngles; ++j) {
        const bool forward = areSameAngle(angles[i], angles[j]);
        const bool reverse = areSameAngle(angles[j], angles[i]);
        EXPECT_EQ(forward, reverse)
            << "symmetry violated for (" << angles[i] << ", " << angles[j] << ")";
    }
}
```

**The wrappers at the bottom**, one line each:

```cpp
TEST(AreSameAngleTest, IdenticalAngles) { asa_test_identical_angles(); }
TEST(AreSameAngleTest, WithinTolerance) { asa_test_within_tolerance(); }
TEST(AreSameAngleTest, OutsideTolerance) { asa_test_outside_tolerance(); }
TEST(AreSameAngleTest, ToleranceBoundary) { asa_test_tolerance_boundary(); }
TEST(AreSameAngleTest, SymmetryAndReflexivity) { asa_test_symmetry_and_reflexivity(); }
TEST(AreSameAngleTest, SignAndWraparoundBehavior) { asa_test_sign_and_wraparound_behavior(); }
TEST(AreSameAngleTest, NonFiniteInputs) { asa_test_non_finite_inputs(); }
```

> One inconsistency to be aware of: this file names its scenarios
> `asa_test_<scenario>()`, but most of the suite uses
> `test_<prefix>_<scenario>()` — about 1170 such functions versus 210. Prefer
> the majority form in new files. Either way, the prefix is what matters.

### The seven scenarios, and why each exists

| Scenario | What it pins down |
| --- | --- |
| `IdenticalAngles` | the ordinary case: a value equals itself |
| `WithinTolerance` | differences below `1E-4` are "same", at small *and* large angles |
| `OutsideTolerance` | differences above `1E-4` are "different" |
| `ToleranceBoundary` | **exactly** `1E-4` is different, because the comparison is strict |
| `SymmetryAndReflexivity` | a property, not a number — survives retuning the tolerance |
| `SignAndWraparoundBehavior` | the *absence* of magnitude compare and `2*pi` reduction |
| `NonFiniteInputs` | NaN and infinity, which the association code can actually supply |

That is the pattern to copy: ordinary case, both sides of every boundary, at
least one property-based check, and explicit coverage of the surprising
behaviour so nobody "fixes" it later.

### Running it

```
$ make unittest
$ ./unit_tests/bin/gtest_main --gtest_filter='AreSameAngleTest.ToleranceBoundary'
[ RUN      ] AreSameAngleTest.ToleranceBoundary

[TEST] asa_test_tolerance_boundary
  Source file:   src/reactions/areSameAngle.cpp
  Function:      areSameAngle(double, double)
  Scenario:      difference is exactly the tolerance 1E-4.
  Pass criteria: strict '<' comparison means false.
    areSameAngle(0, 0.0001) -> false  | expected false  | |diff| = 0.0001  | exactly 1e-4 apart -> not within strict tolerance
    areSameAngle(0.0001, 0) -> false  | expected false  | |diff| = 0.0001  | same boundary case with arguments swapped
    areSameAngle(0, -0.0001) -> false  | expected false  | |diff| = 0.0001  | exactly -1e-4 apart -> also outside
[       OK ] AreSameAngleTest.ToleranceBoundary (0 ms)
[  PASSED  ] 1 test.
```

That console narration is not decoration — it is the house style, and it is
what makes a failure six months from now diagnosable. See §7.

## 5. The GoogleTest you need to know

**`TEST(SuiteName, CaseName) { ... }`** defines a test. gtest finds it
automatically; nothing registers it. In this project the body is always a single
call to a `void test_*()` function defined above it — that keeps the scenarios
as ordinary C++ you can read and reorder.

**`EXPECT_*(...)` checks a condition.** On failure it records the failure and
**keeps going**. The ones you will use:

| You are comparing | Use |
| --- | --- |
| ints, sizes, enums, strings, counters | `EXPECT_EQ(actual, expected)` |
| a double that went through arithmetic | `EXPECT_NEAR(actual, expected, kPfxTol)` |
| a double from a plain assignment | `EXPECT_DOUBLE_EQ(actual, expected)` |
| a bool | `EXPECT_TRUE(x)` / `EXPECT_FALSE(x)` |
| a bound ("stayed in the box") | `EXPECT_LE(x, limit)` / `EXPECT_GE` |

**`ASSERT_*` is the same but aborts the test function immediately.** Default to
`EXPECT_*` so one bad check does not hide the five below it. Use `ASSERT_*` only
when continuing would be *undefined behavior* rather than merely uninformative:

```cpp
// Indexing [1] below would read out of bounds if this failed.
ASSERT_EQ(complexList.size(), 2u) << "dissociation must create a new complex";
```

**Attach a message with `<<`.** Say what the expected behavior *is*, in words —
not a restatement of the numbers, which gtest already prints:

```cpp
EXPECT_EQ(n, 10u) << "table must span bindRadius..Rmax plus one trailing step";
```

**Never use a bare tolerance literal.** Declare one `constexpr double kPfxTol`
at the top of the file and use it throughout, so the precision claim is visible
in one place.

---

## 6. Five rules that keep you from breaking everyone else

All tests share one process and one global namespace. These are not style
preferences.

1. **Prefix every name you define.** Pick a short abbreviation of your source
   file (`szl_`, `rtcs_`, `ailbox_`) and put it on every free function, helper
   and constant. Two files defining `make_molecule()` is a link error for the
   whole suite. Helpers and constants also go inside an unnamed
   `namespace { ... }`.
2. **Do not define `main`, `r`, or `totMatches`.** `gtest_main.cpp` already
   does. Need the RNG? Declare it: `extern gsl_rng* r;`.
3. **Never call anything that can `exit()`, `abort()`, or loop forever.** It
   kills or hangs the entire binary, not just your test. Several NERDSS
   constructors call `exit(1)` on bad input. If a path like that matters, test
   it as a *death test* instead (`EXPECT_EXIT`, with a suite name ending in
   `DeathTest` — see [src/test_error.cpp](src/test_error.cpp)); otherwise note
   it in the file header as not tested and why.
4. **Fully initialize every object you build.** NERDSS routinely indexes into
   vectors without bounds checks and reads members the caller was expected to
   fill in. A half-built `Molecule`, `Complex` or `MolTemplate` segfaults the
   suite. When setup gets long, write one `<prefix>_build_system()` helper and
   call it from every scenario.
5. **Seed the RNG the project's way, never `srand_gsl()`:**

   ```cpp
   if (r == nullptr) {
       const gsl_rng_type* T;
       T = gsl_rng_default;
       r = gsl_rng_alloc(T);
       gsl_rng_set(r, 42);
   }
   ```

   Call that at the top of any scenario that reaches code which might draw a
   random number. Seed 42 always, so results do not change between runs or
   between different `--gtest_filter` selections. Do not assert on a specific
   random draw — assert on an invariant the draw cannot violate.

---

## 7. What to check, and how to narrate it

**Open every scenario with the banner.** Four aligned fields, printed to
`std::cerr`:

```cpp
std::cerr << "\n[TEST] test_szl_basic_table_length\n"
          << "  Source file:   src/reactions/areSameAngle.cpp\n"
          << "  Function:      areSameAngle(double, double)\n"
          << "  Scenario:      <the setup, in one or two lines>\n"
          << "  Pass criteria: <what must hold afterwards>\n";
```

`Pass criteria` is the important one: it is what lets the next person decide
whether a failure is a real regression or a stale test. Then print every input
and result you used to make a decision, the way `asa_check()` prints both
angles, the returned value, the expected value and `|diff|` on one line.

**Cover, for each function in the file:**

- the ordinary case, with a value you worked out by hand from the source;
- each branch you can reach safely;
- the boundary — where a `<=` could plausibly be a `<`;
- degenerate-but-legal inputs the code should survive: empty vectors, empty
  `memberList`, zero-length trajectories, coincident coordinates;
- at least one **invariant** that does not depend on exact numbers ("the complex
  ended up inside the sphere", "a finer step gives a longer table"), so harmless
  numerical changes do not produce false alarms;
- the reversed argument order, if the arguments are symmetric — reactant-ordering
  bugs are a recurring defect class in this codebase.

**Assert against the implementation, not your intuition.** If the code does
something surprising, encode the surprise and explain it in a comment, so nobody
later "fixes" your test. `Coord::operator==` compares components *after*
`roundv()` truncates them to four decimals; a test that picks its own epsilon is
simply wrong about the code.

---

## 8. When it goes wrong

| Symptom | Cause |
| --- | --- |
| `multiple definition of ...` at link time | A name of yours collides with another test file. Prefix it. |
| `undefined reference to <function>` | You included the header but the function lives in a `src/` file that is not compiled, or you misspelled the signature. |
| `No such file or directory` on an include | Use paths relative to `include/`: `"classes/class_Coord.hpp"`. Never `../include/...`, never a bare filename. |
| The whole binary segfaults | Something you constructed is under-initialized (rule 6.4). Run with `--gtest_filter` narrowed to one case to find which. |
| The run hangs | You hit a loop whose step is zero, or a `while` that never advances. |
| Your test passes but checks nothing | Common with file-writing functions: if the output directory does not exist, the `ofstream` fails *silently*. Create the directory, then read the file back and assert on its contents. |

To isolate a crash, narrow the filter until only one case runs:

```sh
./unit_tests/bin/gtest_main --gtest_filter='YourSuite.OneCase'
```

---

## 9. Before you commit

- [ ] `make unittest` builds with no new warnings.
- [ ] `--gtest_filter='YourSuite.*'` passes.
- [ ] The **full** suite still passes — your file is linked into everyone's binary.
- [ ] The file header says what the function does and what you chose not to test.
- [ ] Every name you defined carries your prefix.
- [ ] Any file your test created has been deleted again.

---

## 10. Where to go next

- [STYLE_GUIDE.md](STYLE_GUIDE.md) — the full conventions: naming tables,
  comment density, formatting, death tests, filesystem tests, coverage
  expectations.
- Read before writing: [src/test_class_Coord.cpp](src/test_class_Coord.cpp)
  (a class), [src/test_reflect_traj_check_span.cpp](src/test_reflect_traj_check_span.cpp)
  (a free function with helpers), and
  [src/test_associate_ImplicitLipid_box.cpp](src/test_associate_ImplicitLipid_box.cpp)
  (a large system setup, for when you graduate to association routines).
- [generate_multi_unit_test_gtest.py](generate_multi_unit_test_gtest.py) will
  draft a file for you from a source path. It is a **draft**: it has historically
  needed hand-checking for under-initialized objects that segfault and for
  branches its setup cannot actually reach. Everything above is the acceptance
  criteria for what it gives you.
