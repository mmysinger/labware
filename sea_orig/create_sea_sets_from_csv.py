#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Create original SEA sets from ligand-to-target csv
(Also see sea-build-sets script)

Michael Mysinger 201502 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

import csv
import itertools


class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)
    

def create_sets(infile, sea_root, smiles_column=1, molecules_id_column=2,
                targets_id_column=3, header_lines=1):
    """Create original SEA sets from ligand-to-target csv"""
    field_sep = ";"
    mol_sep = ':'
    in_f = open(infile, "r")
    csv_reader = csv.reader(in_f)
    logging.info("Reading input csv: %s" % infile)
    # Skip header
    for _ in xrange(header_lines):
        csv_reader.next()
    smiles_filename = sea_root + op.extsep + "smi"
    logging.info("Writing SEA smiles file %s" % smiles_filename)
    targets = {}
    with open(smiles_filename, "w") as smiles_f:
        for row in csv_reader:
            smiles = row[smiles_column - 1]
            cid = row[molecules_id_column - 1]
            tid = row[targets_id_column - 1]
            line = smiles + field_sep + cid + "\n"
            smiles_f.write(line)
            clist = targets.setdefault(tid, [])
            clist.append(cid)
    sets_filename = sea_root + op.extsep + "set"
    logging.info("Writing SEA set file %s" % sets_filename)
    with open(sets_filename, "w") as sets_f:
        for k, v in targets.iteritems():
            mols = mol_sep.join(v)
            line = field_sep.join([k, k, mols]) + "\n"
            sets_f.write(line)

def handler(infile, sea_root, **kwargs):
    """I/O handling for the script."""
    try:
        create_sets(infile, sea_root, **kwargs)
    except ScriptError, message:
        logging.error(message)
        return message.value
    return 0

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Create original SEA sets from ligand-to-target csv"
    parser = ArgumentParser(description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("infile", 
                        help="Input ligands to targets map (csv format)")
    parser.add_argument("sea_root", 
                        help="Filename root of SEA output files")
    parser.add_argument("-s", "--smiles-column", type=int, default=1, 
                        help="Column that contains the molecule smiles")
    parser.add_argument("-m", "--molecules-id-column", type=int, default=2, 
                        help="Column that contains the molecule identifiers")
    parser.add_argument("-t", "--targets-id-column", type=int, default=3, 
                        help="Column that contains the target identifiers")
    parser.add_argument("-l", "--header-lines", type=int, default=1, 
                        help="Skip this number of header lines from infile")
    
    options = parser.parse_args(args=argv[1:])
    return handler(options.infile, options.sea_root,
                   smiles_column=options.smiles_column,
                   molecules_id_column=options.molecules_id_column,
                   targets_id_column=options.targets_id_column,
                   header_lines=options.header_lines)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

