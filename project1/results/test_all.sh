#!/bin/bash

# Change into results/ directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

# Make sure all targets are ready
echo "$(. make_all.sh)"

# Clean existing results
rm *.txt

# Run tests
echo "$(. mpitest.sh)"
echo "$(. reductest.sh)"
echo "$(. omptest.sh)"
