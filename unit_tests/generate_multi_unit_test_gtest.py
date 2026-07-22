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
- Use the EXPECT_* assertions provided by googletests and print useful information to stderr.
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

def generate_unit_test(cpp_path: str) -> str:
    with open(cpp_path) as f:
        source = f.read()

    client = anthropic.Anthropic()

    response = client.messages.create(
        model="claude-opus-4-8",
        max_tokens=16384,
        thinking={"type": "adaptive"},
        system=SYSTEM_PROMPT,
        messages=[
            {
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": f"Write a unit test for the file named {cpp_path} with the following content:\n\n```cpp\n{source}\n```",
                        "cache_control": {"type": "ephemeral"},
                    }
                ],
            }
        ],
    )

    return next(block.text for block in response.content if block.type == "text")


if __name__ == "__main__":

    if len(sys.argv) > 1:
        file_name = sys.argv[1]
    else:
        print("You must provide a list of source files to be unit tested.")
        quit()
#    print(f"Generating unit test for: {path}\n")
    with open(file_name, "r") as f:
        for line in f:
            path = "." + line.strip()
            print("Opened file " + path)
            src_path = Path(path)
            out_path = Path("src") / ("test_" + src_path.name)
            out_path.write_text(generate_unit_test(path))
            print(f"Unit test written to: {out_path}")
