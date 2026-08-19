import sys
from pathlib import Path
import anthropic

from dotenv import load_dotenv
load_dotenv()

REPO_ROOT = Path(__file__).resolve().parent.parent
INCLUDE_DIR = REPO_ROOT / "include"
SRC_CLASSES_DIR = REPO_ROOT / "src" / "classes"
TEST_SRC_DIR = REPO_ROOT / "unit_tests" / "src"

# json.hpp is the vendored nlohmann/json single-header (967 KB of the 1.25 MB in
# include/). It contributes nothing to understanding NERDSS APIs, so it is left
# out of the context prefix.
EXCLUDED_HEADERS = {"json.hpp"}

# Existing hand-checked tests used as style references: one class-based, one
# free-function, matching the mix of targets in todo_files.txt.
EXAMPLE_TESTS = ["test_class_Coord.cpp", "test_reflect_traj_check_span.cpp"]

SYSTEM_PROMPT = """\
You are an expert C++ developer. When given a .cpp source file, write a standalone \
unit test for it that includes verbose console output of what function in which source file \
is being tested and what the tests are actually doing.

- Use external test framework gtest.
- Use the EXPECT_* assertions provided by googletests and print useful information to stderr.
- Do not call srand_gsl(); to initialise the random number generator use this instead "const gsl_rng_type *T; T = gsl_rng_default;r = gsl_rng_alloc(T);gsl_rng_set(r, 42);"
- Avoid fatal errors so that all tests run even if some fail.
- Be sure to test all functions in a given file.
- Group related assertions into named void test_*() functions.
- Make sure to add unique prefixes to function names to make sure they don't collide with other test functions.
- Include only the headers needed.
- Make sure the function is commented throughout.
- Have the code be clear about what tests are being run and what criteria is used to pass.
- Do not write a main function because this will be part of a larger suite of tests.
- Do not enclose the source code in a markdown fenced block.
- comment the code of course!
- Return ONLY the C++ source code, no explanation, and do not print any interactive chat.
- Follow the google style guide. """

BUILD_CONVENTIONS = """\
# How these tests are built

Generated tests live in `unit_tests/src/` and are built with `make unittest` from the
repository root. Constraints that follow from that build:

- The compiler standard is **C++17** (`CXXFLAGS = -std=c++17` in the Makefile). Do not use
  language or library features newer than C++17.
- `include/` is on the include path, so project headers are included relative to it:
  `#include "classes/class_Vector.hpp"`, `#include "boundary_conditions/reflect_functions.hpp"`.
  Never use a path relative to the test file (no `../include/...`) and never use a bare
  `#include "class_Vector.hpp"`.
- gtest comes from `pkg-config --cflags/--libs gtest`, so `#include <gtest/gtest.h>` is correct.
- GSL is linked into the test binary, so `<gsl/gsl_rng.h>`, `<gsl/gsl_matrix.h>`,
  `<gsl/gsl_sf_bessel.h>` and friends are available when the code under test needs them.
- `unit_tests/src/gtest_main.cpp` already provides `main()` **and defines these globals**, which
  other translation units declare `extern`:

      gsl_rng* r = nullptr;
      unsigned long totMatches = 0;

  Your file must not define `main`, and must not define `r` or `totMatches` — a duplicate
  definition is a link error for the entire suite. Declare them `extern` if you need them.
- Every generated file is linked into one binary (`unit_tests/bin/gtest_main`), so all
  free function names in your file must be unique across the suite. Prefix them with something
  derived from the file under test.

# NERDSS headers

The complete set of first-party headers from `include/` follows. Use them as the authoritative
source for type definitions, function signatures, default arguments, and struct members — do not
guess at an API that is spelled out below. (The vendored `include/json.hpp` is omitted.)
"""

CLASS_IMPLS_PREAMBLE = """\
# NERDSS class implementations

The definitions for every type declared in `include/classes/` follow. The headers above give you the
signatures; these give you the behavior. Use both to decide what to assert:

- Assert against the behavior the code actually has, including exact rounding and comparison semantics.
  For example, `operator==` on `Coord` does not compare `x`, `y`, `z` directly -- it compares them after
  `roundv`, which truncates to 4 decimal places with a sign-dependent rule. Picking your own epsilon
  instead of reading the implementation produces an assertion that is wrong about the code.
- Do not write a test that reaches an `exit()` or `abort()` path. Several constructors bail out that way
  on bad input (`Coord(std::vector<double>)` calls `exit(1)` when handed more than three values), and
  that kills the entire gtest binary, not just the one case. If a function bails out on bad input, say so
  in a comment rather than triggering it.
- Constructors and member functions frequently index into vectors without bounds checks and read members
  the caller was expected to have filled in first. Read the implementation for those preconditions and
  fully initialize every object your test builds -- an under-initialized `Molecule`, `Complex`,
  `Parameters` or `MolTemplate` segfaults the whole suite.

"""

EXAMPLES_PREAMBLE = """\
# Reference style

Two existing tests from `unit_tests/src/` are shown below. Match their structure, comment
density, and console-output style. They are references for form only — do not copy their
assertions or reuse their test function names.
"""


def collect_headers() -> str:
    """Concatenate every first-party header into one labelled block.

    Sorted deliberately: rglob order is filesystem-dependent, and any reordering
    changes the bytes of the cached prefix and invalidates it.
    """
    files = sorted(INCLUDE_DIR.rglob("*.hpp")) + sorted(INCLUDE_DIR.rglob("*.h"))
    return "\n\n".join(
        f"// ===== include/{h.relative_to(INCLUDE_DIR)} =====\n{h.read_text()}"
        for h in files
        if h.name not in EXCLUDED_HEADERS
    )


def collect_class_impls() -> str:
    """Concatenate every src/classes implementation into one labelled block.

    Sorted for the same reason collect_headers() is: glob order is
    filesystem-dependent, and any reordering changes the bytes of the cached
    prefix and invalidates it.
    """
    return CLASS_IMPLS_PREAMBLE + "\n\n".join(
        f"// ===== src/classes/{c.name} =====\n{c.read_text()}"
        for c in sorted(SRC_CLASSES_DIR.glob("*.cpp"))
    )


def collect_examples() -> str:
    """Concatenate the reference tests into one labelled block."""
    return EXAMPLES_PREAMBLE + "\n\n".join(
        f"// ===== unit_tests/src/{name} =====\n{(TEST_SRC_DIR / name).read_text()}"
        for name in EXAMPLE_TESTS
    )


def build_system_blocks() -> list:
    """Assemble the cached context prefix.

    Ordering is stable -> volatile: these five blocks are byte-identical on every
    request, and the per-file source goes in the user turn after them. Nothing in
    here may vary per file -- no timestamps, no interpolated paths -- or the cache
    never gets read (render order is tools -> system -> messages).

    Two breakpoints, out of the four the API allows. The first one ends the
    expensive headers-plus-implementations run (~130k tokens), so editing
    EXAMPLE_TESTS -- the block most likely to change while tuning the prompt --
    only invalidates the last ~7k tokens. Editing SYSTEM_PROMPT or
    BUILD_CONVENTIONS still invalidates everything, since caching is a prefix
    match and those have to come first.
    """
    return [
        {"type": "text", "text": SYSTEM_PROMPT},
        {"type": "text", "text": BUILD_CONVENTIONS},
        {"type": "text", "text": collect_headers()},
        {"type": "text", "text": collect_class_impls(), "cache_control": {"type": "ephemeral","ttl":"1h"}},
        {"type": "text", "text": collect_examples(), "cache_control": {"type": "ephemeral","ttl":"1h"}},
    ]


def generate_unit_test(client, system_blocks: list, cpp_path: Path) -> str:
    source = cpp_path.read_text()
    rel_path = cpp_path.relative_to(REPO_ROOT)

    # A src/classes target is already in the cached implementations block. It is
    # repeated here anyway -- dropping it per-file would make the prefix vary and
    # kill the cache -- so say why it appears twice.
    note = ("\n\n(This file also appears in the class implementations block above; "
            "it is repeated here as the target.)"
            if cpp_path.parent == SRC_CLASSES_DIR else "")

    # Streamed because thinking is on by default on Opus 5 and max_tokens caps
    # thinking plus output together -- a small non-streaming budget truncates the
    # generated file mid-source.
    count = client.messages.count_tokens(model="claude-opus-5", system=system_blocks,
                             messages=[{"role": "user", "content": "x"}]).input_tokens
    print(f"  input tokens: {count}")
    with client.messages.stream(
        model="claude-opus-5",
        max_tokens=100000,
        thinking={"type": "adaptive"},
        output_config={"effort": "high"},
        system=system_blocks,
        messages=[
            {
                "role": "user",
                "content": f"Write a unit test for {rel_path} with the following content:"
                           f"\n\n```cpp\n{source}\n```{note}",
            }
        ],
    ) as stream:
        response = stream.get_final_message()

    usage = response.usage
    print(
        f"  cache write={usage.cache_creation_input_tokens} "
        f"read={usage.cache_read_input_tokens} "
        f"uncached={usage.input_tokens} out={usage.output_tokens}"
    )

    return next(block.text for block in response.content if block.type == "text")


def read_target_list(list_path: Path) -> list:
    """Read a file list, skipping blank lines and comments.

    todo_files.txt has 22 blank lines; each one used to become the path "."
    and crash the run with IsADirectoryError.
    """
    targets = []
    for line in list_path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        targets.append(REPO_ROOT / line.lstrip("./"))
    return targets


if __name__ == "__main__":
    argv = sys.argv[1:]
    force = "--force" in argv
    positional = [arg for arg in argv if not arg.startswith("-")]

    if not positional:
        print("Usage: generate_multi_unit_test_gtest.py <file-list> [--force]")
        print("  <file-list>  text file listing source paths, e.g. todo_files.txt")
        print("  --force      regenerate tests that already exist")
        quit()

    targets = read_target_list(Path(positional[0]))
    print(f"Read {len(targets)} source files to test.\n")

    print("Assembling repository context (this is cached across files)...")
    system_blocks = build_system_blocks()
    context_chars = sum(len(block["text"]) for block in system_blocks)
    print(f"Context prefix: {context_chars:,} characters\n")

    client = anthropic.Anthropic()

    for src_path in targets:
        out_path = TEST_SRC_DIR / ("test_" + src_path.name)

        if not src_path.is_file():
            print(f"Skipping {src_path}: not found")
            continue
        if out_path.exists() and not force:
            print(f"Skipping {src_path.relative_to(REPO_ROOT)}: "
                  f"{out_path.name} already exists (use --force to regenerate)")
            continue

        print(f"Generating test for {src_path.relative_to(REPO_ROOT)}")
        out_path.write_text(generate_unit_test(client, system_blocks, src_path))
        print(f"  written to {out_path.relative_to(REPO_ROOT)}")
