#!/bin/bash

# zeta1/mach1
mpirun -n 8 ../zeta1/test1
mv vtest.txt zeta1.txt
mpirun -n 8 ../mach1/test1
mv vtest.txt mach1.txt
