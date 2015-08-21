#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Convert SEAware library csv files to SEA original files.

Michael Mysinger 201502 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import csv
import itertools


class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)


def convert_molecules(csv_reader, sea_root):
    """Convert SEAware molecules file to SEA smiles/fingerprints files."""
    logging.info("Reading SEAware molecules file")
    field_sep = ";"
    row = csv_reader.next()
    if row and row[0].startswith("fingerprint_type"):
        logging.info("Skipping fingerprint headers: %s" % row)
        row = csv_reader.next()
    if row and row[0].lower() == "molecule id":
        logging.info("Skipping column headers: %s" % row)
        row = csv_reader.next()
    smiles_filename = sea_root + op.extsep + "smi"
    fingerprints_filename = sea_root + op.extsep + "fp"
    logging.info("Writing SEA smiles file %s" % smiles_filename)
    if len(row) > 2:
        logging.info("Writing SEA fingerprints file %s" %
                     fingerprints_filename)
        fingerprints_found = True
    else:
        fingerprints_found = False
    with open(smiles_filename, "w") as smiles_f, \
          open(fingerprints_filename, "w") as fingerprint_f:
        for row in itertools.chain([row], csv_reader):
            cid, smiles = row[:2]
            if len(row) > 3:
                message =  "Too many fields in molecules file at row: %s" % row
                raise ScriptError(message, 11)
            line = smiles + field_sep + cid + "\n"
            smiles_f.write(line)
            if fingerprints_found:
                fingerprint = row[2]
                line = fingerprint + field_sep + cid + "\n"
                fingerprint_f.write(line)
    if not fingerprints_found:
        os.remove(fingerprints_filename)


def convert_targets(csv_reader, sea_root):
    """Convert SEAware targets file to SEA set file."""
    logging.info("Reading SEAware targets file")
    field_sep = ";"
    row = csv_reader.next()
    if row and row[0].lower() == "target id":
        logging.info("Skipping column headers: %s" % row)
        row = csv_reader.next()
    set_filename = sea_root + op.extsep + "set"        
    logging.info("Writing SEA set file %s" % set_filename)
    with open(set_filename, "w") as set_f:
        for row in itertools.chain([row], csv_reader):
            tid, name, affinity, mols, description = row
            out_tid = name + '_' + affinity
            out_description = description.replace(';', ' -')
            out_description = out_description.replace(':', ' -')
            line = field_sep.join([out_tid, out_description, mols]) + "\n"
            set_f.write(line)
    

def seaware_to_original_sea(molecules, targets, sea_root):
    """Convert SEAware library csv files to SEA original files."""
    molecules_f = open(molecules, "r")
    molecules_reader = csv.reader(molecules_f)
    convert_molecules(molecules_reader, sea_root)

    targets_f = open(targets, "r")
    targets_reader = csv.reader(targets_f)
    convert_targets(targets_reader, sea_root)


def handler(molecules, targets, sea_root):
    """I/O handling for the script."""
    try:
        seaware_to_original_sea(molecules, targets, sea_root)
    except ScriptError, message:
        logging.error(message)
        return message.value
    return 0


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Template: paragraph description goes here."
    parser = ArgumentParser(description=description)
    parser.add_argument("molecules", 
                        help="Input SEAware molecules file (csv format)")
    parser.add_argument("targets", 
                        help="Input SEAware targets file (csv format)")
    parser.add_argument("sea_root", 
                        help="Filename root of SEA output files")
    options = parser.parse_args(args=argv[1:])
    return handler(options.molecules, options.targets, options.sea_root)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

