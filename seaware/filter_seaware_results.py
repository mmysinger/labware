#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Filter SEAware results

Michael Mysinger 201506 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import csv

SELECTED_AFFINITY = "10000"
PVALUE_THRESHOLD = 1.0e-8

module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen
        
def filter_seaware_results(seaware_reader, seaware_writer, best_affinity=False, 
        select_affinity=SELECTED_AFFINITY, pvalue_threshold=PVALUE_THRESHOLD,
        setcore=False):
    """Filter SEAware results by affinity group and p-value."""
    logging.info("Filtering out p-values above: %g" % pvalue_threshold)
    if best_affinity:
        results = {}
        logging.info("Selecting best affinity group")
    else:
        logging.info("Selecting specific affinity group: %s" % select_affinity)
    header = seaware_reader.next()
    logging.info("Skipping input header: %s" % str(header))
    seaware_writer.writerow(header)
    logging.info("Started filtering results")
    selected_counter = 0
    threshold_counter = 0
    for row_counter, row in enumerate(seaware_reader):
        if setcore:
            cid, tuid, affinity, pvalue, maxtc, short, desc = row
        elif len(row) == 8:
            cid, _, tuid, affinity, pvalue, maxtc, short, desc = row
        elif len(row) > 8:
            cid, tuid, affinity, pvalue, maxtc, short, desc = row[:7]
        pvalue = float(pvalue)
        maxtc = float(maxtc)
        if pvalue_threshold and pvalue >= pvalue_threshold:
            threshold_counter += 1
            continue
        if best_affinity:
            key = (cid, tuid)
            if ( key not in results or pvalue < results[key][0] or
                 (maxtc == 1.0 and results[key][1] < 1.0) ):
                results[key] = (pvalue, maxtc, row)
        elif select_affinity and affinity == select_affinity:
            selected_counter += 1
            seaware_writer.writerow(row)
    if best_affinity:
        for key in results:
            selected_counter += 1
            seaware_writer.writerow(results[key][2])
    logging.info("%s result rows read" % (row_counter + 1, ))
    logging.info("%s predictions failed p-value threshold" % threshold_counter)
    logging.info("%s SEAware predictions written" % selected_counter)
    logging.info("Finished filtering results")
        
def handler(in_fn=None, out_fn=None, **kwargs):
    """I/O handling for the script."""
    if in_fn is None:
        in_f = sys.stdin
        in_location = "standard in"
    else:
        in_f = open(in_fn, "r")
        in_location = in_fn
    logging.info("Input SEAware results file: %s" % in_location)
    if out_fn is None:
        out_f = sys.stdout
        out_location = "standard out"
    else:
        out_f = open(out_fn, "w")
        out_location = out_fn
    logging.info("Output filtered results file: %s" % out_location)
    seaware_reader = csv.reader(in_f)
    seaware_writer = csv.writer(out_f)
    try:
        try:
            filter_seaware_results(seaware_reader, seaware_writer, **kwargs)
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
    description = "Filter SEAware results file"
    parser = ArgumentParser(description=description)
    parser.add_argument("-i", "--infile", default=None, 
        help="input SEAware results CSV filename  (default: stdin)")
    parser.add_argument("-o", "--outfile", default=None, 
        help="output filtered results CSV filename (default: stdout)")
    parser.add_argument("-b", "--best-affinity", action="store_true",
                        help="keep only the best affinity set for each target")
    parser.add_argument("-s", "--select-affinity", default=SELECTED_AFFINITY,
                        help="select a particular affinity for every target " +
                        "(default: %(default)s)")
    parser.add_argument("-p", "--pvalue-threshold", default=PVALUE_THRESHOLD,
                        type=float, 
                        help="remove all results worse than this p-value " +
                        "threshold (default: %(default)g)")
    parser.add_argument("--set-mode", action='store_true', 
                        help="Compare set-versus-set results")
    options = parser.parse_args(args=argv[1:])
    return handler(in_fn=options.infile, out_fn=options.outfile,
                   select_affinity=options.select_affinity,
                   best_affinity=options.best_affinity, 
                   pvalue_threshold=options.pvalue_threshold,
                   setcore=options.set_mode)
    options = parser.parse_args(args=argv[1:])
    return handler(infile=options.infile, outfile=options.outfile)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

