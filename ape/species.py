#!/usr/bin/env python3

import logging
import numpy as np

from rmgpy.species import Species, TransitionState

class Species(Species):
    def __init__(self, path=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.path = path

class TransitionState(TransitionState):
    def __init__(self, path=None, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.path = path