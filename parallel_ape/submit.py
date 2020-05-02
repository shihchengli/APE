# -*- coding: utf-8 -*-

"""
Submit scripts
"""
submit_scripts = {
    'parallel_ape': """#!/usr/bin/bash

#PBS -l select=1:ncpus={ncpus}:mpiprocs=1:ompthreads={ncpus}
#PBS -P MST108470
#PBS -q cf40
#PBS -j oe

cd $ape_path/parallel_ape/

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $PBS_JOBID"
echo "=========================================================="

source /home/u1987198/group-sw/bashrc_qchem

export QCSCRATCH=/tmp/ypli/$PBS_JOBID
mkdir -p $QCSCRATCH

source activate ape_env
python Parallel_APE.py {input_file} -n {ncpus} -p {protocol} -mode {sampling_mode} {imaginary_bonds}

rm -r $QCSCRATCH"""}