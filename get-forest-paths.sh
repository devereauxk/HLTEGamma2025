#!/bin/bash

# Check if correct number of arguments provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 FORESTDIR OUTPUT"
    exit 1
fi

FORESTDIR="$1"
OUTPUT="$2"

# Check if FORESTDIR exists
if [ ! -d "$FORESTDIR" ]; then
    echo "Error: Directory $FORESTDIR does not exist"
    exit 1
fi

# Find all .root files recursively and write to output file
find "$FORESTDIR" -name "*.root" -type f > "${OUTPUT}"

echo "Forest paths written to ${OUTPUT}"