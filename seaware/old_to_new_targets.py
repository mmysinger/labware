#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger, SeaChange Pharmaceuticals

Convert SEAware targets file from old (1.7) to new (1.8) format

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

from seacore.util.util import gopen


# Old targets format
#   type_tid, name|chemblid, subtid, cids, description

# New targets format
#   chembl_id, name, subtid, cids, description
        
# Future targets format?
#   tid, subtid, description, cids, user_fields (type, accession, name)


class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

        
def read_accession_map(filename):
    logging.info('\tReading UniProt accession code to ChEMBL id map.')
    accession_map = {}
    accf = gopen(filename)
    accf.next() # skip header
    for line in accf:
        chembl_id, accession = line.split()
        if accession in accession_map:
            mapped_id = accession_map[accession]
            logging.warn("Cannot map accession '%s' " % accession +
                             "to '%s' because it's " % chembl_id +
                             "already mapped to '%s'" % mapped_id)
            continue
        accession_map[accession] = chembl_id
    accf.close()
    return accession_map


def convert_targets(target_reader, target_writer, accession_map):
    """Convert SEAware targets file to new format."""
    logging.info("Converting SEAware targets file")
    field_sep = ";"
    row = target_reader.next()
    if row and row[0].lower() == "target id":
        tid, name, affinity, mols, description = row
        name = "UniProt Name"
        out_row = (tid, name, affinity, mols, description)
        logging.info("Writing new column headers: %s" % (out_row, ))
        target_writer.writerow(out_row)
        row = target_reader.next()
    for row in itertools.chain([row], target_reader):
        tid, name, affinity, mols, description = row
        prefix, suffix = tid.split('_')
        target_type = prefix.upper()
        if target_type == "SP":
            new_id = accession_map[suffix]
        else:  # target_type == "PC" or "PF"
            new_id = name
            name = ""
        out_row = [new_id, name, affinity, mols, description]
        target_writer.writerow(out_row)


def handler(input_targets, output_targets, accession_mapping,
            **kwargs):
    """I/O handling for the script."""
    accession_map = read_accession_map(accession_mapping)
    in_target_f = gopen(input_targets, "rb")
    target_reader = csv.reader(in_target_f)
    out_target_f = gopen(output_targets, "wb")
    target_writer = csv.writer(out_target_f)
    try:
        convert_targets(target_reader, target_writer, accession_map)
    except ScriptError, message:
        logging.error(message)
        return message.value
    finally:
        in_target_f.close()
        out_target_f.close()        
    return 0


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = ( "Convert SEAware targets file from old (1.7) " +
                    "to new (1.8) format" )
    parser = ArgumentParser(description=description)
    parser.add_argument("input_targets", 
                        help="Input SEAware targets file (csv format)")
    parser.add_argument("output_targets", 
                        help="Output SEAware targets file (csv format)")
    parser.add_argument("accession_mapping",
                        help="File mapping ChEMBL Target IDs to " +
                        "Uniprot Accession codes (whitespace separated)")
    options = parser.parse_args(args=argv[1:])
    return handler(options.input_targets, options.output_targets, 
                   options.accession_mapping)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

