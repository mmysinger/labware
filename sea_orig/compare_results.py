#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Compare SEAware results to original SEA (in memory)

Michael Mysinger 201503 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser
import csv
import gzip
import bz2
from itertools import izip_longest

SELECTED_AFFINITY = "10000"
PVALUE_THRESHOLD = 1.0e-5

class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

        
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

def guess_num_comparisons(seaware_value, sea_value):
    """Guess number of comparisons from e-value to p-value ratio."""
    pvalue, _ = seaware_value
    evalue, _ = sea_value
    return evalue/pvalue

def float_close(x, y, threshold=0.5):
    """Check if two floats are within threshold percent of one another."""
    if x > y:
        x, y = y, x
    if (y - x)/y <= threshold/100.0:
        return True
    else:
        return False

def compare_values(seaware_value, sea_value, num_comparisons):
    pvalue, seaware_maxtc = seaware_value
    evalue, sea_maxtc = sea_value
    if float_close(seaware_maxtc, sea_maxtc) and float_close(pvalue * num_comparisons, evalue):
        return True
    else:
        return False

def convert_pvalue(seaware_pvalue, num_comparisons):
    pvalue, maxtc = seaware_pvalue
    return (pvalue * num_comparisons, maxtc)
    
def compare_results(seaware_f, sea_f, num_comparisons=None):
    """Compare SEAware results to original SEA"""
    seaware_f.next() # Skip header
    sea_f.next() # Skip header
    seaware_results = {}
    sea_results = {}
    sea_count = 0
    seaware_count = 0
    yield "Differences between SEAware and SEA files:"
    yield "\t".join(["query_id", "reference_id", "seaware_evalue", "seaware_maxtc", "sea_evalue", "sea_maxtc"]) 
    for seaware_line, sea_line in izip_longest(seaware_f, sea_f):
        if seaware_line:
            cid, smiles, tuid, affinity, pvalue, maxtc, short, desc = seaware_line
            pvalue = float(pvalue)
            maxtc = float(maxtc)
            if not SELECTED_AFFINITY or affinity == SELECTED_AFFINITY:
                if not PVALUE_THRESHOLD or pvalue < PVALUE_THRESHOLD:
                    seaware_count += 1
                    tid = short + "_" + affinity
                    key = (cid, tid)
                    seaware_value = (pvalue, maxtc)
                    if key in sea_results:
                        sea_value = sea_results.pop(key)
                        if num_comparisons is None:
                            num_comparisons = guess_num_comparisons(seaware_value, sea_value) 
                        if not compare_values(seaware_value, sea_value, num_comparisons):
                            yield "\t".join(str(x) for x in key + convert_pvalue(seaware_value, num_comparisons) + sea_value)
                    else:
                        seaware_results[key] = seaware_value
        if sea_line:
            cid, tid, evalue, maxtc = sea_line
            evalue = float(evalue)
            maxtc = float(maxtc)
            if not SELECTED_AFFINITY or tid.endswith("_" + SELECTED_AFFINITY):
                if not PVALUE_THRESHOLD or not num_comparisons or evalue < (num_comparisons * PVALUE_THRESHOLD):
                    sea_count += 1
                    key = (cid, tid)
                    sea_value = (evalue, maxtc)
                    if key in seaware_results:
                        seaware_value = seaware_results.pop(key)
                        if num_comparisons is None:
                            num_comparisons = guess_num_comparisons(seaware_value, sea_value) 
                        if not compare_values(seaware_value, sea_value, num_comparisons):
                            yield "\t".join(str(x) for x in key + convert_pvalue(seaware_value, num_comparisons) + sea_value)
                    else:
                        sea_results[key] = sea_value
    # Clean up SEA results passed threshold before num_comparisons was guessed
    for key, value in sea_results.items():
        evalue, maxtc = value
        if PVALUE_THRESHOLD and num_comparisons and evalue >= (num_comparisons * PVALUE_THRESHOLD):
            sea_results.pop(key)
            sea_count -= 1
    yield ""
    yield "Read %d seaware and %d sea results that passed thresholds." % (seaware_count, sea_count)
    yield "Used %.1f set comparisons to convert p-values to e-values." % num_comparisons
    if len(seaware_results):
        yield ""
        yield "Missing matching results for the following seaware results:"
        for key, value in seaware_results.iteritems():
            yield "\t".join(str(x) for x in key + convert_pvalue(value, num_comparisons))
    if len(sea_results):
        yield ""
        yield "Missing matching results for the following sea results:"
        for key, value in sea_results.iteritems():
            yield "\t".join(str(x) for x in key + value)

def handler(seaware_fn, sea_fn, out_fn=None, **kwargs):
    """I/O handling for the script."""
    if out_fn is None:
        out_f = sys.stdout
    else:
        out_f = gopen(out_fn, "w")
    out_f.write("Using SEAware file: %s\n" % seaware_fn)
    seaware_f = gopen(seaware_fn , "r")
    seaware_reader = csv.reader(seaware_f)
    out_f.write("Using SEA file: %s\n" % sea_fn)
    sea_f = gopen(sea_fn , "r")
    sea_reader = csv.reader(sea_f)
    
    try:
        try:
            for result in compare_results(seaware_reader, sea_reader, **kwargs):
                out_f.write(result + "\n")
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        seaware_f.close()
        sea_f.close()
        out_f.close()
    return 0


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Compare SEAware results to original SEA."
    parser = ArgumentParser(description=description)
    parser.add_argument("seaware", 
                        help="Input SEAware results file (csv format)")
    parser.add_argument("sea", 
                        help="Input original SEA scores file (csv format)")
    parser.add_argument("-n", "--num-comparisons", type=int, 
                        help="Number of set comparisons printed during original SEA run (default: approximate guess)")
    parser.add_argument("-o", "--outfile",
                        help="Output filename (default: stdout)")
    options = parser.parse_args(args=argv[1:])
    return handler(options.seaware, options.sea, out_fn=options.outfile, num_comparisons=options.num_comparisons)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

