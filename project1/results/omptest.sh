#!/bin/bash

# Change into results/ directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"


# zeta2/mach2

# 2 Threads
OMP_NUM_THREADS=2
../zeta2/test2
mv vtest.txt zeta2_t2.txt
../mach2/test2
mv vtest.txt mach2_t2.txt

# 4 Threads
OMP_NUM_THREADS=4
../zeta2/test2
mv vtest.txt zeta2_t4.txt
../mach2/test2
mv vtest.txt mach2_t4.txt

# 8 Threads
OMP_NUM_THREADS=8
../zeta2/test2
mv vtest.txt zeta2_t8.txt
../mach2/test2
mv vtest.txt mach2_t8.txt

# 16 Threads
OMP_NUM_THREADS=16
../zeta2/test2
mv vtest.txt zeta2_t16.txt
../mach2/test2
mv vtest.txt mach2_t16.txt
