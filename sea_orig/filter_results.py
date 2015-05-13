#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Filter original SEA results

Michael Mysinger 201503 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser
import csv

SELECTED_AFFINITY = "10000"
EVALUE_THRESHOLD = 1.0e-3

module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen

def filter_results(in_csv, out_csv, select_affinity=SELECTED_AFFINITY,
                   best_affinity=False, evalue_threshold=EVALUE_THRESHOLD):
    """Filter original SEA results"""
    header = in_csv.next()
    logging.info("Skipping input header: %s" % str(header))
    out_csv.writerow(header)
    if best_affinity:
        results = {}
    logging.info("Filtering results")
    for row in in_csv:
        query_id, reference_id, evalue, maxtc = row
        evalue = float(evalue)
        if evalue_threshold and evalue >= evalue_threshold:
            continue
        try:
            target_id, affinity = reference_id.rsplit("_", 1)
        except ValueError:
            target_id = reference_id
            affinity = select_affinity
        if best_affinity:
            key = (query_id, target_id)
            if key not in results or results[key][0] > evalue:
                results[key] = (evalue, row)
        elif select_affinity and select_affinity == affinity:
            out_csv.writerow(row)
    if best_affinity:
        logging.info("Selecting best affinity result")
        for key in results:
            out_csv.writerow(results[key][1])
    logging.info("Finished")

def handler(in_fn=None, out_fn=None, **kwargs):
    """I/O handling for the script."""
    if in_fn is None:
        in_f = sys.stdin
        logging.info("Inputing scores file from standard in")
    else:
        in_f = gopen(in_fn, "r")
        logging.info("Inputing scores file from %s" % in_fn)
    if out_fn is None:
        out_f = sys.stdout
        logging.info("Outputing scores file to standard out")
    else:
        out_f = gopen(out_fn, "w")
        logging.info("Outputing scores file to %s" % out_fn)
    in_csv = csv.reader(in_f)
    out_csv = csv.writer(out_f)
    try:
        try:
            filter_results(in_csv, out_csv, **kwargs)
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
    description = "Compare SEAware results to original SEA."
    parser = ArgumentParser(description=description)
    parser.add_argument("-i", "--infile", default=None, 
                        help="input SEA scores filename (default: stdin)")
    parser.add_argument("-o", "--outfile", default=None, 
                        help="output filename (default: stdout)")
    parser.add_argument("-b", "--best-affinity", action="store_true",
                        help="keep only the best affinity set for each target")
    parser.add_argument("-s", "--select-affinity", default=SELECTED_AFFINITY,
                        help="select a particular affinity for every target " +
                        "(default: %(default)s)")
    parser.add_argument("-e", "--evalue-threshold", default=EVALUE_THRESHOLD,
                        type=float, 
                        help="remove all results worse than this e-value " +
                        "threshold (default: %(default)g)")
    options = parser.parse_args(args=argv[1:])
    return handler(in_fn=options.infile, out_fn=options.outfile,
                   select_affinity=options.select_affinity,
                   best_affinity=options.best_affinity,
                   evalue_threshold=options.evalue_threshold)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

