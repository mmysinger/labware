#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Convert SEA original files to SEAware library csv files.

Michael Mysinger 201502 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import csv
import itertools
from collections import OrderedDict

class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)


def convert_molecules(sea_root, csv_writer):
    """Convert SEA smiles/fingerprints files to SEAware molecules file."""
    field_sep = ";"
    smiles_filename = sea_root + op.extsep + "smi"
    logging.info("Reading SEA smiles file %s" % smiles_filename)
    molecules = OrderedDict()
    with open(smiles_filename, 'r') as smiles_f:
        for line in smiles_f:
            smiles, cid = line.strip().split(field_sep)
            if cid in molecules:
                logging.error("Molecule id %s with smiles %s overwritten" %
                              (cid, molecules[cid]))
            molecules[cid] = (smiles, None)

    fingerprints_filename = sea_root + op.extsep + "fp"
    logging.info("Reading SEA fingerprints file %s" % fingerprints_filename)
    with open(fingerprints_filename, 'r') as fingerprints_f:
        for line in fingerprints_f:
            fp, cid = line.strip().split(field_sep)
            if cid in molecules:
                smiles, _ = molecules[cid] 
                molecules[cid] = (smiles, fp)
            else:
                logging.error("Molecule id %s with fp %s has no smiles" %
                              (cid, fp))

    logging.info("Writing SEAware molecules file")
    csv_writer.writerow(["molecule id", "smiles", "fingerprint"])
    for cid, (smiles, fp) in molecules.iteritems():
        csv_writer.writerow([cid, smiles, fp])


def convert_targets(sea_root, csv_writer):
    """Convert SEA set file to SEAware targets file."""
    field_sep = ";"
    set_filename = sea_root + op.extsep + "set"        
    logging.info("Reading SEA set file %s" % set_filename)
    logging.info("Writing SEAware targets file")
    csv_writer.writerow(["target id", "name", "affinity", "molecules",
                         "description"])
    with open(set_filename, 'r') as set_f:
        for line in set_f:
            target, description, mols = line.strip().split(field_sep)
            try:
                tid, affinity = target.rsplit("_", 1)
            except ValueError:
                tid = target
                affinity = ""
            name = tid
            csv_writer.writerow([tid, name, affinity, mols, description])
    

def sea_to_seaware(molecules, targets, sea_root):
    """Convert SEA original files to SEAware library csv files."""
    
    molecules_f = open(molecules, "w")
    molecules_writer = csv.writer(molecules_f)
    convert_molecules(sea_root, molecules_writer)

    targets_f = open(targets, "w")
    targets_writer = csv.writer(targets_f)
    convert_targets(sea_root, targets_writer)


def handler(molecules, targets, sea_root):
    """I/O handling for the script."""
    try:
        sea_to_seaware(molecules, targets, sea_root)
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
    parser.add_argument("sea_root", 
                        help="Filename root of SEA input files")
    parser.add_argument("molecules", 
                        help="Output SEAware molecules file (csv format)")
    parser.add_argument("targets", 
                        help="Output SEAware targets file (csv format)")
    options = parser.parse_args(args=argv[1:])
    return handler(options.molecules, options.targets, options.sea_root)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

