#!/bin/bash

#PBS -N serial_n512
#PBS -q training
#PBS -W group_list=imf_lille-tma4280
#PBS -l walltime=00:05:00
#PBS -l select=1:ncpus=20

cd $PBS_O_WORKDIR
module load gcc
module load openmpi
module load openblas

exec ../../../build/serial 512
