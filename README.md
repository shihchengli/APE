# Automated Property Estimator | APE

[![MIT license](http://img.shields.io/badge/license-MIT-brightgreen.svg)](http://opensource.org/licenses/MIT)

**Automated Property Estimator (APE)** software is a tool to automatically calculate the thermodynamic properties with considering anharmonic effect. All calculations were done using a development version of the [QChem](http://www.q-chem.com/) software package.

Currently, APE can only estimate the thermodynamic property of a given molecule.

**This software is still work in progress to develop more functions.**

## Clone and setup path

Befor installing APE, download [Anacinda](anaconda.com/download/) and [Git](https://git-scm.com/downloads).

Install the latest version of APE by cloning the source code via Git.

- `git clone https://github.com/shihchengli/APE.git`

Add APE to your local path in .bashrc (make sure to change ~/Path/to/APE/ accordingly):

- `export PYTHONPATH=$PYTHONPATH:~/Path/to/APE/`

## Install dependencies

- Install the latest DEVELOPER version of RMG (which has Arkane). It is recommended to follow RMG’s [Developer installation by source using Anaconda](http://reactionmechanismgenerator.github.io/RMG-Py/users/rmg/installation/index.html#for-developers-installation-by-source-using-anaconda-environment) instructions. Make sure to add RMG-Py to your PATH and PYTHONPATH variables as explained in RMG’s documentation. 

- Clone [ARC](https://reactionmechanismgenerator.github.io/ARC/installation.html) and [pysisyphus](https://github.com/eljost/pysisyphus) in a separate folder and add it to your PYTHONPATH just as well.

  `git clone https://github.com/ReactionMechanismGenerator/ARC.git`
  
  `export PYTHONPATH=$PYTHONPATH:~/Path/to/ARC/`

  `git clone https://github.com/eljost/pysisyphus.git`
  
  `export PYTHONPATH=$PYTHONPATH:~/Path/to/pysisyphus/`

- Create an Anaconda environment for APE:

  `conda env create -f environment.yml`
  
- Finally either close and reopen your terminal to refresh your environment variables, or type the following command.

  `source ~/.bashrc`

- Activate the APE environment every time before you run APE:

  `source activate ape_env`

## How to run a job

To run APE, make sure to first activate the APE environment.

- Basic run

  `python <path_to_the_APE_folder>/APE.py freq_output.q.out -n <ncpus> -p <UMN or UMVT (default)> > thermo_result`

- Using the API
  ```markdown

  from ape.main import APE
  
  ape = APE(input_file, name=None, project_directory=None, protocol=None, multiplicity=None, charge = None, external_symmetry=None, level_of_theory=None, basis=None, ncpus=None)
  
  ape.execute()
  ```
## How to give feedback
  
If you are having issues, please post to our <a href="https://github.com/shihchengli/APE/issues">Issues</a> page. We will do our best to assist.

## Credits

- [Professor Yi-Pei Li's research group](https://webpageprodvm.ntu.edu.tw/Li-group/Default.aspx) at [National Taiwan University](https://www.ntu.edu.tw/english/spotlight/2020/1820_20200327.html)
