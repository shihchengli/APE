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


    group = parser.add_mutually_exclusive_group()
    group.add_argument('-p', type=str, help='the sampling protocol(UMN or UMVT) chossen (default: UMVT)')

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

    ape_object = APE(input_file = input_file, project_directory=project_directory, protocol=protocol)
    ape_object.execute()

################################################################################

if __name__ == '__main__':
    main()
