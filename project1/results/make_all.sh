#!/bin/bash

# Change into results/ directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"

# Run make in all test-relevant directories
cd '../zeta1'
make
cd '../zeta2'
make
cd '../mach1'
make
cd '../mach2'
make
cd '../reduc'
make
