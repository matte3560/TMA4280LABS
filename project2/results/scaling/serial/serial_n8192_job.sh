#!/bin/bash

#PBS -N serial_n8192
#PBS -q training
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=20

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
module load openblas

exec ../../../build/serial 8192
