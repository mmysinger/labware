#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger, SeaChange Pharmaceuticals

Convert SEAware molecules file from old (1.7) to new (1.8) format

Michael Mysinger 201508 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import bz2
import csv
import gzip
import itertools

from fpcore import inchi_key
from seacore.util.util import gopen


# Old molecules format
#   cid, smiles, fp

# New molecules format
#   cid, smiles, fp, user fields (inchi_key)


class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)


def convert_molecules(mol_reader, mol_writer):
    """Convert SEAware molecules file to new format."""
    logging.info("Converting SEAware molecules file")
    field_sep = ";"
    row = mol_reader.next()
    if row and row[0].startswith("fingerprint_type"):
        logging.info("Copying fingerprint headers: %s" % row)
        mol_writer.writerow(row)
        row = mol_reader.next()
    if row and row[0].lower() == "molecule id":
        fp_header = "Fingerprint"
        if len(row) > 2 and row[2]:
            fp_header = row[2]
        out_row = [row[0], row[1], fp_header, "InChi Key"]
        logging.info("Writing new column headers: %s" % out_row)
        mol_writer.writerow(out_row)
        row = mol_reader.next()
    if len(row) > 2:
        logging.info("Copying fingerprints column")
        fingerprints_found = True
    else:
        fingerprints_found = False
    for row in itertools.chain([row], mol_reader):
        cid, smiles = row[:2]
        if len(row) > 3:
            message =  "Too many fields in molecules file at row: %s" % row
            raise ScriptError(message, 11)
        ikey = inchi_key(smiles)
        if not ikey:
            logging.warn("Skipping molecule with SMILES '%s'" % smiles)
            continue
        fingerprint = ""
        if fingerprints_found:
            fingerprint = row[2]
        mol_writer.writerow([cid, smiles, fingerprint, ikey])        


def handler(input_molecules, output_molecules, **kwargs):
    """I/O handling for the script."""
    in_mol_f = gopen(input_molecules, "rb")
    mol_reader = csv.reader(in_mol_f)
    out_mol_f = gopen(output_molecules, "wb")
    mol_writer = csv.writer(out_mol_f)
    try:
        convert_molecules(mol_reader, mol_writer)
    except ScriptError, message:
        logging.error(message)
        return message.value
    finally:
        in_mol_f.close()
        out_mol_f.close()
    return 0


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = ( "Convert SEAware molecules file from old (1.7) " +
                    "to new (1.8) format" )
    parser = ArgumentParser(description=description)
    parser.add_argument("input_molecules", 
                        help="Input SEAware molecules file (csv format)")
    parser.add_argument("output_molecules", 
                        help="Output SEAware molecules file (csv format)")
    options = parser.parse_args(args=argv[1:])
    return handler(options.input_molecules, options.output_molecules)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

