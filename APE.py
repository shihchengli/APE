# -*- coding: utf-8 -*-

"""
APE - Automated Property Estimator
"""

import argparse
import os

from ape.main import APE

def parse_command_line_arguments(command_line_args=None):
    
    parser = argparse.ArgumentParser(description='Automated Property Estimator (APE)')
    parser.add_argument('file', metavar='FILE', type=str, nargs=1,
                        help='a file describing the job to execute')
    parser.add_argument('-p', type=str, help='the sampling protocol(UMN or UMVT) chossen (default: UMVT)')
    parser.add_argument('-i', type=str, help='the imaginary bonds for QMMM calculation')

    args = parser.parse_args(command_line_args)
    args.file = args.file[0]

    return args

def main():
    """ The main APE executable function"""
    args = parse_command_line_arguments()
    input_file = args.file
    protocol = args.p
    project_directory = os.path.abspath(os.path.dirname(args.file))
    if not protocol:
        protocol = 'UMVT'
        print('This calculation will use UMVT as sampling protocol')
    elif protocol == 'UMN' or protocol == 'UMVT':
        print('This calculation will use {} as sampling protocol'.format(protocol))
    
    # imaginary bonds for QMMM calculation
    # atom indices starts from 1
    imaginary_bonds = args.i
    if args.i is not None:
        imaginary_bonds_string = imaginary_bonds.strip('[').strip(']')
        imaginary_bonds = []
        for bond in imaginary_bonds_string.split(','):
            atom1, atom2 = bond.split('-')
            imaginary_bonds.append([int(atom1), int(atom2)])

    ape_object = APE(input_file = input_file, project_directory=project_directory, protocol=protocol, imaginary_bonds=imaginary_bonds)
    ape_object.execute()

################################################################################

if __name__ == '__main__':
    main()
