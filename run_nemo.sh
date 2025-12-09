#!/bin/bash --login

#PBS -N nemo
#PBS -l walltime=01:00:00
#PBS -q collab
#PBS -l select=1
#PBS -P other

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

ulimit -c unlimited
ulimit -s unlimited

NEMO_N=64
XIOS_N=8

echo " mpiexec --cpu-bind=depth -n $NEMO_N -d 1 -n $XIOS_N -d 1 ./nemo"
mpiexec --cpu-bind=depth -n $NEMO_N -d 1 -n $XIOS_N -d 1 ./nemo
