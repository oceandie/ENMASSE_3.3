#!/bin/bash --login

#PBS -N domcfg
#PBS -l walltime=00:30:00
#PBS -q collab
#PBS -l select=1
#PBS -P other

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

ulimit -c unlimited
ulimit -s unlimited

NEMO_N=4
XIOS_N=1

echo " mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./make_domain_cfg.exe"
mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./make_domain_cfg.exe
