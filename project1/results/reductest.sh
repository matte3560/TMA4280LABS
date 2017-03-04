#!/bin/bash

# reduc size1
mpirun -n 2 ../reduc/zeta 8388608
mv result.txt reduc_n2_s1.txt
mpirun -n 4 ../reduc/zeta 8388608
mv result.txt reduc_n4_s1.txt
mpirun -n 8 ../reduc/zeta 8388608
mv result.txt reduc_n8_s1.txt
# reduc_rd size1
mpirun -n 2 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_n2_s1.txt
mpirun -n 4 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_n4_s1.txt
mpirun -n 8 ../reduc/zeta_rd 8388608
mv result.txt reduc_rd_n8_s1.txt

# reduc size2
mpirun -n 2 ../reduc/zeta 16777216
mv result.txt reduc_n2_s2.txt
mpirun -n 4 ../reduc/zeta 16777216
mv result.txt reduc_n4_s2.txt
mpirun -n 8 ../reduc/zeta 16777216
mv result.txt reduc_n8_s2.txt
# reduc_rd size2
mpirun -n 2 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_n2_s2.txt
mpirun -n 4 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_n4_s2.txt
mpirun -n 8 ../reduc/zeta_rd 16777216
mv result.txt reduc_rd_n8_s2.txt
