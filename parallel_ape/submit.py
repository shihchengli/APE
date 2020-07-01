# -*- coding: utf-8 -*-

"""
Submit scripts
"""
submit_scripts = {
    'parallel_ape': """#!/usr/bin/bash

#PBS -l select=1:ncpus={ncpus}:mpiprocs=1:ompthreads={ncpus}
#PBS -q workq
#PBS -j oe


cd {job_path}

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $PBS_JOBID"
echo "=========================================================="

module load qchem

export QCSCRATCH=/tmp/$PBS_JOBID
mkdir -p $QCSCRATCH

source ~/tools/anaconda3/etc/profile.d/conda.sh
conda activate ape_env
python $ape_path/parallel_ape/Parallel_APE.py {input_file} -n {ncpus} -p {protocol} -mode {sampling_mode} {imaginary_bonds}

rm -r $QCSCRATCH"""}