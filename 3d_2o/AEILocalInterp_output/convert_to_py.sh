#!/bin/bash

# ====================================================
# Script Name: convert.sh
# Description: Converts all [base].compute.c files in a specified directory
#              to trusted_[base].py files by applying specific transformations
#              and prepending a Python header.
# Usage: ./convert.sh /path/to/directory
# ====================================================

# Exit immediately if a command exits with a non-zero status
set -e

# Function to display usage information
usage() {
    echo "Usage: $0 /path/to/directory"
    exit 1
}

# Check if exactly one argument (the directory path) is provided
if [ "$#" -ne 1 ]; then
    echo "Error: Exactly one argument expected."
    usage
fi

# Assign the first argument to DIRECTORY variable
DIRECTORY="$1"

# Check if the provided argument is a valid directory
if [ ! -d "$DIRECTORY" ]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# Define the Python header to prepend to each .py file
PYTHON_HEADER='import sympy as sp


def RATIONAL(x, y):
    return sp.Rational(x, y)


x, y, z = sp.symbols("x y z")
'

# Enable nullglob to handle cases where no files match the pattern
shopt -s nullglob

# Initialize a counter for processed files
count=0

# Loop over all [base].compute.c files in the specified directory
for c_file in "$DIRECTORY"/*.compute.c; do
    # Since nullglob is set, the loop does not iterate if no matches are found

    # Increment the counter
    count=$((count + 1))

    # Extract the base filename without '.compute.c' (e.g., coeffs_dxx.compute.c -> coeffs_dxx)
    base_name=$(basename "$c_file" .compute.c)

    # Replace hyphens with underscores in base_name
    base_name_underscored=${base_name//-/_}

    # Define the output .py filename (e.g., trusted_coeffs_dxx.py)
    py_file="$DIRECTORY/trusted_${base_name_underscored}.py"

    echo "Processing '$c_file' -> '$py_file'..."

    # Create or overwrite the .py file with the Python header
    echo "$PYTHON_HEADER" > "$py_file"

    # Apply the transformations:
    # 1. Exclude lines starting with 'fp '.
    # 2. Remove leading six spaces.
    # 3. Join lines until a semicolon is found at the end.
    # 4. Replace '->' with '__'.
    # 5. Remove semicolons.
    # 6. Replace 'coeffs__' with 'trustedcoeffs__'.
    #
    # These transformations are chained using a single sed command for efficiency.
    grep -v '^fp ' "$c_file" | \
    sed 's/^      //g; :a; /;$/!{ N; s/\n/ /; ba }; s/->/__/g; s/;//g; s/coeffs_/trusted_coeffs_/g' >> "$py_file"

    # Inform the user of successful processing
    echo "Successfully processed '$c_file' -> '$py_file'."
done

# Disable nullglob after processing
shopt -u nullglob

# Check if any files were processed
if [ "$count" -eq 0 ]; then
    echo "No *.compute.c files found in directory '$DIRECTORY'."
else
    echo "Successfully processed $count file(s)."
fi
