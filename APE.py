# -*- coding: utf-8 -*-

"""
APE - Automated Property Estimator
"""

import os
import logging

from ape.main import APE

ape = APE()

# Parse and validate the command-line arguments
args = ape.parse_command_line_arguments()

# Execute the job
ape.execute()

try:
    import psutil

    process = psutil.Process(os.getpid())
    memory_info = process.memory_info()
    logging.info('Memory used: %.2f MB' % (memory_info.rss / 1024.0 / 1024.0))
except ImportError:
    logging.info('Optional package dependency "psutil" not found; memory profiling information will not be saved.')