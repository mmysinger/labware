#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Template: one line summary goes here.

Michael Mysinger 201505 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

# XXX - MMM currently set for bitterdb, may want better defaults 
CUTOFF_PVALUE = 1.0e-7 # 1e-4 original
CUTOFF_EF = 2.0        # Nat2012: EF > 1
CUTOFF_MINPAIRS = 3 # Nat2012: Target-ADR pairs > 10 retained

class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)
        
# reverse a dictionary of dictionaries
def flip_dict(x2drugs):
    drugs2x = {}
    for x, drugs in x2drugs.iteritems():
        for drug in drugs:
            drugs2x[drug] = drugs2x.get(drug, {})
            drugs2x[drug][x] = True
    return drugs2x

def read_target_predictions(seaware_reader):
    pass

def ef_analysis(in_f):
    """Module code starts here."""
    yield "Hello World!\n"

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
            for result in ef_analysis(in_f, **kwargs):
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
    description = "Compute enrichment factors."
    parser = ArgumentParser(description=description)
    parser.add_argument("seaware_results",  
                        help="SEAware results")
    parser.add_argument("-o", "--outfile", default=None, 
                        help="output file (default: stdout)")
    options = parser.parse_args(args=argv[1:])
    return handler(infile=options.infile, outfile=options.outfile)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

