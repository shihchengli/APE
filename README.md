# Automated Property Estimator | APE

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

**Automated Property Estimator (APE)** software is a tool to automatically calculate the thermodynamic properties with considering anharmonic effect. All calculations were done using a development version of the [QChem](http://www.q-chem.com/) software package.

Currently, APE can estimate the thermodynamic property of a given molecule and the kinetics properties of a given elementary reaction in gas-phase.

Please check out our published paper of [Comparative Analysis of Uncoupled Mode Approximations for Molecular Thermochemistry and Kinetics](https://doi.org/10.1021/acs.jctc.2c00664) in the Journal of Chemical Theory and Computation.

## Clone and setup path

Befor installing APE, download [Anacinda](anaconda.com/download/) and [Git](https://git-scm.com/downloads).

Install the latest version of APE by cloning the source code via Git.

- `git clone https://github.com/shihchengli/APE.git`

Add APE to your local path in .bashrc (make sure to change ~/Path/to/APE/ accordingly):

- `export PYTHONPATH=$PYTHONPATH:~/Path/to/APE/`

## Install dependencies

- Install the v3.1.0 version of RMG (which has Arkane). It is recommended to follow RMG’s [Developer installation by source using Anaconda](http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/index.html#for-developers-installation-by-source-using-anaconda-environment) instructions. Make sure to add RMG-Py to your PATH and PYTHONPATH variables as explained in RMG’s documentation.

- Clone the forked repo of [ARC](https://reactionmechanismgenerator.github.io/ARC/installation.html) from my github page in a separate folder and add it to your PYTHONPATH just as well.

  `git clone https://github.com/shihchengli/ARC.git`
  
  `export PYTHONPATH=$PYTHONPATH:~/Path/to/ARC/`

- Create an Anaconda environment for APE:

  `conda env create -f environment.yml`
  
- Finally either close and reopen your terminal to refresh your environment variables, or type the following command.

  `source ~/.bashrc`

- Activate the APE environment every time before you run APE:

  `source activate ape_env`

## Add APE aliases to your .bashrc (for convenience)

  ```markdown

  export ape_path=$HOME'/Path/to/APE/'
  alias ape_env='source activate ape_env'
  alias ape_dir='cd $ape_path'

  ```

## How to run a job

To run APE, make sure to first activate the APE environment.

- Basic run
  
  `python <path_to_the_APE_folder>/APE.py -n <ncpus> input.py`

- Using the API
  ```markdown

  from ape.main import APE
  
  ape = APE(input_file, project_directory=None)
  
  ape.execute()
  ```
## How to give feedback
  
If you are having issues, please post to our <a href="https://github.com/shihchengli/APE/issues">Issues</a> page. We will do our best to assist.

## Licence notifications

This project contains source material from the geomeTRIC package.

Copyright 2016-2020 Regents of the University of California and the Authors

Authors: Lee-Ping Wang, Chenchen Song\ https://github.com/leeping/geomeTRIC/blob/master/LICENSE

## How to cite
Please include the following citation if APE was used for an academic study:

- Shih-Cheng Li, Yen-Chun Lin, Yi-Pei Li, Comparative Analysis of Uncoupled Mode Approximations for Molecular Thermochemistry and Kinetics, J. Chem. Theory Comput. 2022, XXXX, XXX, XXX-XXX (DOI: [10.1021/acs.jctc.2c00664](https://doi.org/10.1021/acs.jctc.2c00664))

## Credits

- [Professor Yi-Pei Li's research group](https://webpageprodvm.ntu.edu.tw/Li-group/Default.aspx) at [National Taiwan University](https://www.ntu.edu.tw/english/spotlight/2020/1820_20200327.html)
