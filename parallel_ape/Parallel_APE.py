# -*- coding: utf-8 -*-

"""
Parallel_APE
"""

import argparse
import os

from parallel_ape.parallel_main import Parallel_APE

def parse_command_line_arguments(command_line_args=None):
    
    parser = argparse.ArgumentParser(description='Automated Property Estimator (APE)')
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
                        help='a file describing the job to execute')
    parser.add_argument('-n', type=int, help='number of CPUs to run quantum calculation')
    parser.add_argument('-p', type=str, default='UMVT', help='the sampling protocol(UMN or UMVT) chossen (default: UMVT)')
    parser.add_argument('-i', type=str, help='the imaginary bonds for QMMM calculation')
    parser.add_argument('-mode', type=int, help='the specific mode to sample points')

    args = parser.parse_args(command_line_args)
    args.file = args.file[0]

    return args

def main():
    """ The main parallel APE executable function"""
    args = parse_command_line_arguments()
    input_file = args.file
    ncpus = args.n
    protocol = args.p
    project_directory = os.path.abspath(os.path.dirname(args.file))
    mode = args.mode
    # imaginary bonds for QMMM calculation
    # atom indices starts from 1
    imaginary_bonds = args.i
    if args.i is not None:
        imaginary_bonds_string = imaginary_bonds.strip('[').strip(']')
        imaginary_bonds = []
        for bond in imaginary_bonds_string.split(','):
            atom1, atom2 = bond.split('-')
            imaginary_bonds.append([int(atom1), int(atom2)])


    ape_object = Parallel_APE(input_file = input_file, ncpus=ncpus, project_directory=project_directory, protocol=protocol, imaginary_bonds=imaginary_bonds)
    ape_object.parse()
    ape_object.sampling(sampling_mode=mode)

################################################################################

if __name__ == '__main__':
    main()
