#!/bin/bash

#PBS -N mpi_p16_n512
#PBS -q training
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=20:mpiprocs=16

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
module load openblas

mpirun ../../../build/mpi 512
