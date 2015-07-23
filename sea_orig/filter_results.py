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

def filter_results(sea_reader, sea_writer, select_affinity=SELECTED_AFFINITY,
                   best_affinity=False, evalue_threshold=EVALUE_THRESHOLD):
    """Filter original SEA results"""
    header = sea_reader.next()
    logging.info("Skipping input header: %s" % str(header))
    sea_writer.writerow(header)
    if best_affinity:
        results = {}
    logging.info("Filtering results")
    for row in sea_reader:
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
            if ( key not in results or results[key][0] > evalue or
                 (float(maxtc) == 1.0  and float(results[key][1][3]) < 1.0) ):
                results[key] = (evalue, row)
        elif select_affinity and select_affinity == affinity:
            sea_writer.writerow(row)
    if best_affinity:
        logging.info("Selecting best affinity results")
        for key in results:
            sea_writer.writerow(results[key][1])
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
    sea_reader = csv.reader(in_f)
    sea_writer = csv.writer(out_f)
    try:
        try:
            filter_results(sea_reader, sea_writer, **kwargs)
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
    description = "Filter original SEA results"
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

