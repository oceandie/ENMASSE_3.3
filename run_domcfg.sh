#!/bin/bash --login

#PBS -N domcfg
#PBS -l walltime=00:10:00
#PBS -q normal
#PBS -l select=1:ncpus=256:coretype=genoa
#PBS -P other

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

ulimit -c unlimited
ulimit -s unlimited

NEMO_N=1
XIOS_N=1

echo " mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./make_domain_cfg.exe"
mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./make_domain_cfg.exe
