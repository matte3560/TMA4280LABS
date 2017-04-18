#!/bin/bash

#PBS -N hybrid_p4t9
#PBS -q training
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:05:00
#PBS -l select=2:ncpus=20:mpiprocs=2:ompthreads=9

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
module load openblas

mpirun ../../build/hybrid 16384
