#!/bin/bash --login

#PBS -N nemo
#PBS -l walltime=01:00:00
#PBS -q collab
#PBS -l select=1:coretype=milan
#PBS -P other

source $DATADIR/NEMO/ukmo_utils/load_hpc_modules_xios3.sh

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

ulimit -c unlimited
ulimit -s unlimited

NEMO_N=30
XIOS_N=2

echo " mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./nemo -n $XIOS_N -d 1 ./xios_server.exe"
mpiexec --cpu-bind=depth -n $NEMO_N -d 1 ./nemo -n $XIOS_N -d 1 ./xios_server.exe
