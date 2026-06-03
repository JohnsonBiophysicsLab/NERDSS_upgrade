import sys
import anthropic

from dotenv import load_dotenv
load_dotenv()

SYSTEM_PROMPT = """\
You are an expert C++ developer. When given a .cpp source file, write a standalone \
unit test for it in the same style as this example test:

- No external test framework (no gtest, no catch2) — plain C++ with a main().
- Use helper functions like require_close(actual, expected, label) and \
require_true(condition, label) that print to stderr and call std::exit(1) on failure.
- Group related assertions into named void test_*() functions.
- Include only the headers needed.
- Return ONLY the C++ source code, no explanation."""

def generate_unit_test(cpp_path: str) -> str:
    with open(cpp_path) as f:
        source = f.read()

    client = anthropic.Anthropic()

    response = client.messages.create(
        model="claude-opus-4-8",
        max_tokens=2048,
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
    path = sys.argv[1] if len(sys.argv) > 1 else "./src/reactions/norm_function.cpp"
    print(f"Generating unit test for: {path}\n")
    print(generate_unit_test(path))
