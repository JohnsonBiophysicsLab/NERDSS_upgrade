import anthropic

from dotenv import load_dotenv
load_dotenv()  # add before anthropic.Anthropic()

client = anthropic.Anthropic()
messages=[
        {"role": "user", "content": "What is 2 + 2?"}
    ]
response = client.messages.create(
    model="claude-opus-4-8",
    max_tokens=256,
    messages=messages
)
print(messages[0]["content"])
print("here's my response!")
for block in response.content:
    if block.type == "text":
        print(block.text)
