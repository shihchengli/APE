# -*- coding: utf-8 -*-
"""
This module contains classes which extend Exception for usage in the APE module.
"""

class ConvergeError(Exception):
    """
    An exception raised when the convergence is not met.
    """
    pass

class InputError(Exception):
    pass

class JobError(Exception):
    """
    An exception class for exceptional behavior that occurs while working with jobs.
    """
    pass