#!/usr/bin/env python

import math

# User inputs -------------------------
jpni = 4
jpnj = 5
nemo_nprocpernode = 32
xios_nproc = 8
xios_nprocpernode = 2
# -------------------------------------

hpc_queue = "normal"
nemo_nproc = jpni * jpnj # we don't use land suppression 
nemo_nodes = math.ceil(nemo_nproc / nemo_nprocpernode)
xios_nodes = math.ceil(xios_nproc / xios_nprocpernode)
tot_nodes = math.ceil(nemo_nodes + xios_nodes)

with open ('run_nemo.sh', 'w') as rsh:
    rsh.write(f"""\
#! /bin/bash --login

#PBS -N nemo
#PBS -l walltime=02:00:00
#PBS -q {hpc_queue}
#PBS -l select={tot_nodes}:ncpus=256:coretype=genoa
#PBS -P other

source ~/NEMO/load_hpc_modules_xios3.sh

export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)
export OMP_NUM_THREADS=1
cd $PBS_O_WORKDIR

ulimit -c unlimited
ulimit -s unlimited

echo "mpiexec -n {nemo_nproc} -ppn {nemo_nprocpernode} --cpu-bind=list:$(./bind_list.py {nemo_nprocpernode} 192) ./nemo : -n {xios_nproc} -ppn {xios_nprocpernode} --cpu-bind=list:$(./bind_list.py {xios_nprocpernode} 192) xios_server.exe"

mpiexec -n {nemo_nproc} -ppn {nemo_nprocpernode} --cpu-bind=list:$(./bind_list.py {nemo_nprocpernode} 192) ./nemo : -n {xios_nproc} -ppn {xios_nprocpernode} --cpu-bind=list:$(./bind_list.py {xios_nprocpernode} 192) xios_server.exe
""")
