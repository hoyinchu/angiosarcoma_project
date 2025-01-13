#!/bin/bash

# Ensure the script receives the correct number of arguments
if [ "$#" -ne 3 ]; then
  echo "Usage: $0 <input_file> <keyword_file> <output_file>"
  exit 1
fi

input_file="$1"
keyword_file="$2"
output_file="$3"

# Check if keyword_file exists
if [ ! -f "$keyword_file" ]; then
  echo "Keyword file '$keyword_file' not found!"
  exit 1
fi

# Create a grep pattern from the keyword file
patterns=$(awk '{print $0}' "$keyword_file" | paste -sd '|' -)

# Use grep with the pattern to filter lines containing any of the keywords
grep -E "$patterns" "$input_file" > "$output_file"

echo "Lines containing any keywords from '$keyword_file' have been written to '$output_file'."
