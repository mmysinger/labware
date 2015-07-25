#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Create random event clusters for enrichment analysis

Michael Mysinger 201506 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import csv
import itertools
import random

import numpy as np
from sklearn.neighbors import KernelDensity
from sklearn.grid_search import GridSearchCV

module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen


DEFAULT_NUM_CLUSTERS = 100


def silvermans_rule(data, verbose=False):
    sstddev = np.std(data, ddof=1)
    bandwidth_estimate = (4.0/3/len(data))**0.2 * sstddev
    if verbose:
        mean = np.mean(data)
        logging.info("Mean: %g" % mean)
        logging.info("Sample standard deviation: %g" % sstddev)
        logging.info("Silverman's bandwidth estimate: %g" % bandwidth_estimate)
    return bandwidth_estimate


def select_kde(data, crossfold=20):
    bw_est = silvermans_rule(data, verbose=True)
    logging.info("Selecting bandwidth by %d crossfold validation" % crossfold)
    grid = GridSearchCV(KernelDensity(),
        {"bandwidth": np.linspace(0.1 * bw_est, 2.0 * bw_est, 39)},
        cv=crossfold)
    grid.fit(data[:, None])
    kde = grid.best_estimator_
    logging.info("Bandwidth selected by crossfold validation: %g" %
                 kde.bandwidth) 
    return kde


def read_target_sizes(targets_reader):
    field_sep = ";"
    mols_sep = ":"
    logging.info("Reading SEAware targets file")
    row = targets_reader.next()
    if row and row[0] == "target id":
        logging.info("Skipping column headers: %s" % row)
        row = targets_reader.next()
    for row in itertools.chain([row], targets_reader):
        tid, name, affinity, mols, description = row
        yield mols.count(mols_sep) + 1


def read_molecule_ids(in_reader):
    logging.info("Reading SEAware input CSV file")
    field_sep = ";"
    row = in_reader.next()
    if row and row[0].startswith("fingerprint_type"):
        logging.info("Skipping fingerprint headers: %s" % row)
        row = in_reader.next()
    if row and row[0] == "molecule id":
        logging.info("Skipping column headers: %s" % row)
        row = in_reader.next()
    cids = set()
    for row in itertools.chain([row], in_reader):
        cid, smiles = row[:2]
        if len(row) > 3:
            message =  "Too many fields in molecules file at row: %s" % row
            raise ScriptError(message, 11)
        cids.add(cid)
    return cids


def seaware_input_to_random_events(in_reader, targets_reader,
                                   num_clusters=DEFAULT_NUM_CLUSTERS):
    """Create random event clusters for enrichment analysis"""
    np.random.seed(42)
    random.seed(42)
    sizes = list(read_target_sizes(targets_reader))
    kde = select_kde(np.array(sizes))
    cids = read_molecule_ids(in_reader)
    for i in xrange(num_clusters):
        random_size = int(round(kde.sample()[0][0]))
        while random_size < 1:
            random_size = int(round(kde.sample()[0][0]))
        print random_size
        cluster_id = "cluster_%04d" % (i+1)
        for cid in random.sample(cids, random_size):
            yield [cid, cluster_id, "", cluster_id]


def handler(in_fn, targets_fn, out_fn, **kwargs):
    """I/O handling for the script."""
    logging.info("Input file: %s" % in_fn)
    in_f = gopen(in_fn)
    in_reader = csv.reader(in_f)
    logging.info("Targets file: %s" % targets_fn)
    targets_f = gopen(targets_fn)
    targets_reader = csv.reader(targets_f)
    logging.info("Output file: %s" % out_fn)
    out_f = gopen(out_fn, "wb")
    out_writer = csv.writer(out_f)
    try:
        try:
            for result in seaware_input_to_random_events(
                    in_reader, targets_reader, **kwargs):
                out_writer.writerow(result)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        in_f.close()
        targets_f.close()
        out_f.close()
    return 0

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Create random event clusters for enrichment factor analysis. The SEAware 'input' CSV file contains a collection of negative molecules that will serve as the random background. The SEAware 'targets' CSV file is used to estabilish the reference distribution of target sizes. Random samples are then drawn using a non-parametric gaussian kernel density estimator, where the bandwidth is choosen using crossfold validation. The output is an events file containing the desired number of random event clusters."
    parser = ArgumentParser(description=description)
    parser.add_argument("input",
                        help="SEAware input CSV file")
    parser.add_argument("targets",
                        help="SEAware targets CSV file")
    parser.add_argument("output", 
                        help="output events CSV file")
    parser.add_argument("-n", "--num-clusters", default=DEFAULT_NUM_CLUSTERS,
        help="number of random events clusters (default: %(default)s)")
    options = parser.parse_args(args=argv[1:])
    return handler(in_fn=options.input, targets_fn=options.targets,
                   out_fn=options.output, num_clusters=options.num_clusters)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

