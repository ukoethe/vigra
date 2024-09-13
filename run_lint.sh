#!/bin/bash

# Find all .c, .h, and .cpp files in the current directory
files=$(find . -name '*.c' -o -name '*.h' -o -name '*.cpp')

# Initialize a flag to track if any trailing whitespace is found
error_found=0

# Loop through each file and check for trailing whitespace
for file in $files; do
    if grep -q "[[:blank:]]$" "$file"; then
        echo "Trailing whitespace found in: $file"
        error_found=1
    fi
done

# If trailing whitespace was found, exit with an error
if [ $error_found -eq 1 ]; then
    exit 1
else
    echo "No trailing whitespace found."
fi
