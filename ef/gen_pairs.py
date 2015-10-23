#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Compute pairs of targets for target-target to event EF analysis.

Michael Mysinger 201510 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import csv
from collections import defaultdict

module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen


CUTOFF_MINPAIRS = 4       # Nat2012: Target-ADR pairs > 10 retained

Target = namedtuple("Target", "name description")


def read_results(results_reader):
    """Read targets to molecules mapping"""
    logging.info("Reading targets")
    header = results_reader.next()
    logging.info("Skipping SEAware results header: %s" % str(header))
    molecules = {}
    targets = {}
    targets_to_drugs = defaultdict(set)
    for row in results_reader:
        if len(row) == 8:
            cid, smiles, tid, affinity, pvalue, maxtc, name, desc = row
            mol_extra = []
        elif len(row) > 8:
            cid, tid, affinity, pvalue, maxtc, name, desc, smiles = row[:8]
            mol_extra = row[8:]
        if cid not in molecules:
            molecules[cid] = (smiles, mol_extra)
        if tid not in targets:
            targets[tid] = Target(name, desc)
        # implicitly takes the union over remaining affinity groups
        targets_to_drugs[tid].add(cid, affinity, pvalue, maxtc)
    logging.info("Read %d targets" % len(targets_to_drugs))
    return targets_to_drugs, molecules, targets


def gen_pairs(results_reader, min_pairs=CUTOFF_MINPAIRS):
    """Compute pairs of targets for target-target to event EF analysis"""
    logging.info("Using min-pairs cutoff = %d" % min_pairs)
    targets_to_drugs, molecules, targets = read_results(results_reader)

# Psuedo-code
#   Create target -> cid mapping with sea results along for the ride
#   Consider all target pairs (sort pair name alphabetically?)
#       If number of common drugs in target pair > min_pairs
#           Write out those drugs results




def handler(events_fn, results_fn, out_fn, **kwargs):
    """I/O handling for the script."""
    logging.info("SEAware results file: %s" % results_fn)
    results_f = open(results_fn, "r")
    results_reader = csv.reader(results_f)
    out_f = open(out_fn, "w")
    logging.info("Output file: %s" % out_fn)
    out_writer = csv.writer(out_f)
    try:
        try:
            for result in gen_pairs(results_reader, **kwargs):
                out_writer.writerow(result)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        results_f.close()
        out_f.close()
    return 0


def main(argv):
    """Parse arguments."""
    log_level = logging.INFO
    log_format = "%(levelname)s: %(message)s"
    logging.basicConfig(level=log_level, format=log_format)
    description = "Compute pairs of targets for target-target to event EF analysis."
    parser = ArgumentParser(description=description)
    parser.add_argument("results",  
                        help="SEAware results mapping molecules to targets")
    parser.add_argument("output", 
                        help="output target pairs CSV file")
    parser.add_argument("-m", "--min-pairs", type=int, default=CUTOFF_MINPAIRS, 
        help="Minimum pairs cutoff for target pair creation (default: %(default)s)")
    options = parser.parse_args(args=argv[1:])
    # Add file logger
    log_fn = options.output.replace(".csv", "") + ".log"
    file_handler = logging.FileHandler(log_fn, mode="w")
    log_formatter = logging.Formatter(log_format)
    file_handler.setFormatter(log_formatter)
    file_handler.setLevel(log_level)
    root_logger = logging.getLogger()
    root_logger.addHandler(file_handler)
    return handler(results_fn=options.results, 
                   out_fn=options.output, min_pairs=options.min_pairs)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

