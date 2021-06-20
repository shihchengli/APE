# -*- coding: utf-8 -*-

"""
APE's main module.
To run APE through its API, first get a freq output file from Qchem, then call the .execute. For example:

    ape = APE(input_file)
    ape.execute()

"""

import argparse
import logging
import os
import sys
import time

from ape.input import load_input_file
from ape.sampling import SamplingJob
from ape.thermo import ThermoJob
from ape.kinetics import KineticsJob

class APE(object):
    """
    The main APE class.
    """

    def __init__(self, input_file=None, output_directory=None):
        self.job_list = []
        self.input_file = input_file
        self.output_directory = output_directory
    
    def parse_command_line_arguments(self):
        parser = argparse.ArgumentParser(description='Automated Property Estimator (APE)')
        parser.add_argument('file', metavar='FILE', type=str, nargs=1, 
                            help='a file describing the job to execute')
        parser.add_argument('-n', '--ncpus',type=int, help='number of CPUs to run quantum calculation')

        # Options for controlling the amount of information printed to the console 
        parser.add_argument('-q', '--quiet', action='store_const', const=logging.WARNING, default=logging.INFO,
                           dest='verbose', help='only print warnings and errors')
        parser.add_argument('-v', '--verbose', action='store_const', const=logging.DEBUG, default=logging.INFO,
                           dest='verbose', help='print more verbose output')
        parser.add_argument('-d', '--debug', action='store_const', const=0, default=logging.INFO, dest='verbose',
                           help='print debug information')

        # Add options for whether running this calculation in parallel or not
        parser.add_argument('-parallel', '--run-in-parallel', type=bool, default=False, help='run parallel jobs (default: False)')

        # Add options for controlling what directories files are written to
        parser.add_argument('-o', '--output-directory', type=str, nargs=1, default='',
                            metavar='DIR', help='use DIR as output directory')

        # Add options for controlling generation of plots
        parser.add_argument('-p', '--plot',type=bool, default=True, help='generating plots')

        args = parser.parse_args()

        # Extract the input file
        self.input_file = args.file[0]

        # Extract the number of cpus
        self.ncpus = args.ncpus

        # Extract the log verbosity
        self.verbose = args.verbose

        # Extract the APE running version (serial or parallel) 
        self.run_in_parallel = args.run_in_parallel

        # Extract the plot settings
        self.plot = args.plot

        # Determine the output directory
        # By default the directory containing the input file is used, unless an
        # alternate directory is specified using the -o flag
        if args.output_directory and os.path.isdir(args.output_directory[0]):
            self.output_directory = os.path.abspath(args.output_directory[0])
        else:
            self.output_directory = os.path.dirname(os.path.abspath(args.file[0]))
        return args

    def load_input_file(self, input_file):
        """
        Load a set of jobs from the given `input_file` on disk. Returns the
        loaded set of jobs as a list.
        """
        self.input_file = input_file
        self.job_list, self.reaction_dict, self.species_dict, \
            self.transition_state_dict = load_input_file(self.input_file, self.output_directory)
        return self.job_list

    def execute(self):
        """
        Execute, in order, the jobs found in input file specified by the
        `input_file` attribute.
        """
        # Initialize the logging system (both to the console and to a file in the
        # output directory)
        initialize_log(self.verbose, os.path.join(self.output_directory, 'ape.log'))

        # Print some information to the beginning of the log
        log_header(run_in_parallel=self.run_in_parallel, ncpus=self.ncpus)

        # Load the input file for the job
        self.job_list = self.load_input_file(self.input_file)

        # Initialize (and clear!) the output files for the job
        if self.output_directory is None:
            self.output_directory = os.path.dirname(os.path.abspath(self.input_file))
        output_file = os.path.join(self.output_directory, 'output.py')
        with open(output_file, 'w'):
            pass

        # run SamplingJob
        for job in self.job_list:
            if isinstance(job, SamplingJob):
                job.ncpus = self.ncpus
                job.execute()

        # run thermo and kinetics jobs
        for job in self.job_list:
            if isinstance(job, ThermoJob):
                job.ncpus = self.ncpus
                job.load_save()
                job.execute()
            if isinstance(job, KineticsJob):
                job.execute(output_directory=self.output_directory, plot=self.plot)

        # Print some information to the end of the log
        log_footer()

def initialize_log(verbose=logging.INFO, log_file=None):
    """
    Set up a logger for APE to use to print output to stdout. The
    `verbose` parameter is an integer specifying the amount of log text seen
    at the console; the levels correspond to those of the :data:`logging` module.
    """
    # Create logger
    logger = logging.getLogger()
    logger.setLevel(verbose)

    # Use custom level names for cleaner log output
    logging.addLevelName(logging.CRITICAL, 'Critical: ')
    logging.addLevelName(logging.ERROR, 'Error: ')
    logging.addLevelName(logging.WARNING, 'Warning: ')
    logging.addLevelName(logging.INFO, '')
    logging.addLevelName(logging.DEBUG, '')
    logging.addLevelName(0, '')

    # Create formatter and add to handlers
    formatter = logging.Formatter('%(levelname)s%(message)s')

    # Remove old handlers before adding ours
    while logger.handlers:
        logger.removeHandler(logger.handlers[0])

    # Create console handler; send everything to stdout rather than stderr
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(verbose)
    ch.setFormatter(formatter)
    logger.addHandler(ch)

    # Create file handler; always be at least verbose in the file
    if log_file:
        fh = logging.FileHandler(filename=log_file)
        fh.setLevel(min(logging.DEBUG, verbose))
        fh.setFormatter(formatter)
        logger.addHandler(fh)


def log_header(level=logging.INFO, run_in_parallel=None, ncpus=None):
    """
    Output a header containing identifying information about APE to the log.
    """
    logging.log(level, 'APE execution initiated at {0}'.format(time.asctime()))
    if run_in_parallel and ncpus is not None:
        logging.log(level, 'Parallel version, running on {0} processors'.format(ncpus))
    elif ncpus is not None:
        logging.log(level, 'Serial version, running on {0} processors'.format(ncpus))
    logging.log(level, '')
    logging.log(level, '################################################################')
    logging.log(level, '#                                                              #')
    logging.log(level, '#              Automated Property Estimator (APE)              #')
    logging.log(level, '#                                                              #')
    logging.log(level, '#   Developer: Shih-Cheng Li (f08524007@ntu.edu.tw)            #')
    logging.log(level, '#   P.I.:      Yi-Pei Li (yipeili@ntu.edu.tw)                  #')
    logging.log(level, '#   Website:   https://webpageprodvm.ntu.edu.tw/Li-group/      #')
    logging.log(level, '#                                                              #')
    logging.log(level, '################################################################')
    logging.log(level, '')


def log_footer(level=logging.INFO):
    """
    Output a footer to the log.
    """
    logging.log(level, '')
    logging.log(level, 'APE execution terminated at {0}'.format(time.asctime()))