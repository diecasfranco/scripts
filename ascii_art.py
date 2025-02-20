#!/usr/bin/env python3

import pyfiglet

# Get user input
input_text = input("Enter your text: ")
input_emoji = input("Enter your emoji: ")

# Generate ASCII art
ascii_art = pyfiglet.figlet_format(input_text)

# Add emoji to ASCII art
if input_emoji:
    emoji_art = f"{input_emoji} " * len(input_text)
    ascii_art = f"{ascii_art}\n{emoji_art}"

# Print ASCII art
print(ascii_art)
