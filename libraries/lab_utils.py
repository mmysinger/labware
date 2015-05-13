#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Miscellaneous reusable library functions

Michael Mysinger 201504 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

def libraries_path():
    """Example for how to bootstrap this package into a script"""
    module_path = os.path.realpath(os.path.dirname(__file__)) 
    labware_path = os.path.join(module_path, "..")
    sys.path.append(labware_path)
    #from libraries.lab_utils import ScriptError, gopen


def gopen(fname, mode='rU'):
    """Flexibly open compressed files."""
    if os.path.splitext(fname)[-1] == '.%s' % 'gz':
        if 'b' not in mode:
            # gzip does not support universal newlines
            mode = mode.replace('U', '') + 'b'
        f = gzip.open(fname, mode)
    elif os.path.splitext(fname)[-1] == '.%s' % 'bz2':
        if 'b' not in mode:
            # bzip2 does support universal newlines
            mode = mode + 'b'
        f = bz2.BZ2File(fname, mode)
    else:
        f = open(fname, mode)
    return f


def float_close(x, y, threshold=0.5):
    """Check if two floats are within threshold percent of one another."""
    if x > y:
        x, y = y, x
    if (y - x)/y <= threshold/100.0:
        return True
    else:
        return False
