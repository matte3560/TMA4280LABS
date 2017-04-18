#!/bin/bash

#PBS -N vtest_p1t1
#PBS -q training
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=20:mpiprocs=1:ompthreads=1

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
module load openblas

mpirun ../../../build/hybrid_vtest
