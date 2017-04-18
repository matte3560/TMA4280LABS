#!/bin/bash

#PBS -N vtest_p36t4
#PBS -q training
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:05:00
#PBS -l select=2:ncpus=20:mpiprocs=18:ompthreads=4

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
module load openblas

mpirun ../../../build/hybrid_vtest
