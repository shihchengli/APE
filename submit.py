# -*- coding: utf-8 -*-

"""
Submit scripts
"""

#!/usr/bin/bash

#PBS -l select=1:ncpus=8:mpiprocs=1:ompthreads=8
#PBS -P MST108470
#PBS -q cf40
#PBS -j oe

submit_scripts = {
    'ape': """#!/usr/bin/bash

#PBS -l select=1:ncpus={cpus}:mpiprocs=1:ompthreads={cpus}
#PBS -P MST108470
#PBS -q cf40
#PBS -j oe

cd /home/r08524007/qchem

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $PBS_JOBID"
echo "=========================================================="

source /home/u1987198/group-sw/bashrc_qchem

export QCSCRATCH=/tmp/ypli/$PBS_JOBID
mkdir -p $QCSCRATCH

python APE.py -nt {cpus}  {input} {output}

rm -r $QCSCRATCH"""