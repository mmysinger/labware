#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Generate fingerprints for original SEA using seaware-academic FPCore

Michael Mysinger 201503 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

from fpcore.rdkit_common import clean_smiles
from fpcore.fingerprint import fingerprint_g, FINGERPRINT_TYPES

class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

def load_smiles(in_f):
    field_sep = ';'
    for line in in_f:
        smi, cid = line.strip().split(field_sep)
        yield cid, smi, None

def generate_fingerprints(in_f, fp_type=FINGERPRINT_TYPES[0]):
    """Generate fingerprints for original SEA using seaware-academic FPCore"""
    field_sep = ';'
    clean_smiles_g = clean_smiles(load_smiles(in_f))
    fp_g = fingerprint_g(((smi, cid) for cid, smi, _ in clean_smiles_g),
                         fp_type=fp_type) 
    for fp, cid in fp_g:
        yield fp + field_sep + cid + '\n'


def handler(infile=None, outfile=None, **kwargs):
    """I/O handling for the script."""
    if infile is None:
        in_f = sys.stdin
    else:
        in_f = open(infile, "r")
    if outfile is None:
        out_f = sys.stdout
    else:
        out_f = open(outfile, "w")
    try:
        try:
            for result in generate_fingerprints(in_f, **kwargs):
                out_f.write(result)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        in_f.close()
        out_f.close()
    return 0

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Generate fingerprints for original SEA using seaware-academic FPCore"
    parser = ArgumentParser(description=description)
    parser.add_argument("-i", "--infile", default=None, 
                        help="input file (default: stdin)")
    parser.add_argument("-o", "--outfile", default=None, 
                        help="output file (default: stdout)")
    parser.add_argument("-f", "--fingerprint", default=FINGERPRINT_TYPES[0], 
                        help="fingerprint type (default: %(default)s)")
    options = parser.parse_args(args=argv[1:])
    return handler(infile=options.infile, outfile=options.outfile,
                   fp_type=options.fingerprint)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

