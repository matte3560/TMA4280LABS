#!/bin/bash

# Change into results/ directory
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$DIR"


# reduc size1
mpirun -n 2 ../reduc/zeta 8388608
mv result.txt reduc_p2_s1.txt
mpirun -n 4 ../reduc/zeta 8388608
mv result.txt reduc_p4_s1.txt
mpirun -n 8 ../reduc/zeta 8388608
mv result.txt reduc_p8_s1.txt
mpirun -n 16 ../reduc/zeta 8388608
mv result.txt reduc_p16_s1.txt
mpirun -n 32 ../reduc/zeta 8388608
mv result.txt reduc_p32_s1.txt

# reduc_rd size1
mpirun -n 2 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_p2_s1.txt
mpirun -n 4 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_p4_s1.txt
mpirun -n 8 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_p8_s1.txt
mpirun -n 16 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_p16_s1.txt
mpirun -n 32 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_p32_s1.txt


# reduc size2
mpirun -n 2 ../reduc/zeta 16777216
mv result.txt reduc_p2_s2.txt
mpirun -n 4 ../reduc/zeta 16777216
mv result.txt reduc_p4_s2.txt
mpirun -n 8 ../reduc/zeta 16777216
mv result.txt reduc_p8_s2.txt
mpirun -n 16 ../reduc/zeta 16777216
mv result.txt reduc_p16_s2.txt
mpirun -n 32 ../reduc/zeta 16777216
mv result.txt reduc_p32_s2.txt

# reduc_rd size2
mpirun -n 2 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_p2_s2.txt
mpirun -n 4 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_p4_s2.txt
mpirun -n 8 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_p8_s2.txt
mpirun -n 16 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_p16_s2.txt
mpirun -n 32 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_p32_s2.txt
