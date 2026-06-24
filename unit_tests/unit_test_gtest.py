import sys
from pathlib import Path
import anthropic

from dotenv import load_dotenv
load_dotenv()

SYSTEM_PROMPT = """\
You are an expert C++ developer. When given a .cpp source file, write a standalone \
unit test for it that includes verbose console output of what function in which source file \
is being tested and what the tests are actually doing.

- Use external test framework gtest.
- Use helper functions like require_close(actual, expected, label) and \
require_true(condition, label) that print to stderr and call std::exit(1) on failure.
- Group related assertions into named void test_*() functions.
- Include only the headers needed.
- Make sure the function is commented throughout.
- Have the main code be clear about what tests are being run and what criteria is used to pass.
- comment the code of course!
- Return ONLY the C++ source code, no explanation."""

def generate_unit_test(cpp_path: str) -> str:
    with open(cpp_path) as f:
        source = f.read()

    client = anthropic.Anthropic()

    response = client.messages.create(
        model="claude-opus-4-8",
        max_tokens=4096,
        system=SYSTEM_PROMPT,
        messages=[
            {
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": f"Write a unit test for this file:\n\n```cpp\n{source}\n```",
                        "cache_control": {"type": "ephemeral"},
                    }
                ],
            }
        ],
    )

    return next(block.text for block in response.content if block.type == "text")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        print("You must provide a source file to be unit tested.")
        quit()
    print(f"Generating unit test for: {path}\n")

    src_path = Path(path)
    out_path = Path("src") / ("test_" + src_path.name)

    out_path.write_text(generate_unit_test(path))
    print(f"Unit test written to: {out_path}")
