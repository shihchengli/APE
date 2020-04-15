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
    """
    An exception raised when doing QMMM calculation without specifing imaginary bonds to describe frustrated translation and rotation.
    """
    pass