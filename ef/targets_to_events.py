#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Create events mapping from SEAware targets csv file for EF analysis

Michael Mysinger 201506 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import csv
import itertools


module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen

        
def targets_to_events(targets_reader):
    """Module code starts here."""
    field_sep = ";"
    mols_sep = ":"
    tid_sep = "_"
    logging.info("Reading SEAware targets file")
    row = targets_reader.next()
    if row and row[0] == "target id":
        logging.info("Skipping column headers: %s" % row)
        row = targets_reader.next()
    for row in itertools.chain([row], targets_reader):
        tid, name, affinity, mols, description = row
        out_id = tid
        if affinity:
            out_id = tid + tid_sep + affinity 
        for mol in mols.split(mols_sep):
            yield [mol, out_id]
    logging.info("Finished writing events file")


def handler(in_fn=None, out_fn=None, **kwargs):
    """I/O handling for the script."""
    if in_fn is None:
        in_f = sys.stdin
        in_location = "standard in"
    else:
        in_f = gopen(in_fn)
        in_location = in_fn
    logging.info("Input file: %s" % in_location)
    if out_fn is None:
        out_f = sys.stdout
        out_location = "standard out"
    else:
        out_f = gopen(out_fn, "w")
        out_location = out_fn
    logging.info("Output file: %s" % out_location)
    targets_reader = csv.reader(in_f)
    events_writer = csv.writer(out_f)
    try:
        try:
            for result in targets_to_events(targets_reader, **kwargs):
                events_writer.writerow(result)
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
    description = ( "Create events mapping from SEAware targets csv file " +
                    "for EF analysis" )
    parser = ArgumentParser(description=description)
    parser.add_argument("-i", "--infile", default=None, 
                        help="input file (default: stdin)")
    parser.add_argument("-o", "--outfile", default=None, 
                        help="output file (default: stdout)")
    options = parser.parse_args(args=argv[1:])
    return handler(in_fn=options.infile, out_fn=options.outfile)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

