#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Generate a heatmap from EF q-values CSV file

Adapted from original code by Leo Gendelev

Michael Mysinger 201504 Created
Michael Mysinger 201508 Modified for EF analysis
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser
import csv
import numpy
import matplotlib
from matplotlib import pyplot
import seaborn

module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen
from sea_orig.clustering import hierarchical_x, hierarchical_y, hierarchical_both

QVALUE_CUTOFF = 0.001
EF_CUTOFF = 6.0


def build_array(value_dict, x_labels, y_labels):
    values = []
    for y_label in y_labels:
        values.append([value_dict.get((x_label, y_label), 0.0)
                            for x_label in x_labels])
    return numpy.array(values)


def plot_heatmap(xy_data, out_fn, x_labels=None, y_labels=None):
    logging.info("Input data shape: %s" % str(xy_data.shape))
    figure, axes = pyplot.subplots()
    x = max(7 + 0.15 * len(xy_data[0]), 10)
    y = max(5 + 0.15 * len(xy_data), 10)
    logging.info("Using figure size: %.1f x %.1f" %  (x, y))
    figure.set_size_inches(x, y)
    if x_labels == None:
        x_labels = ""
    if y_labels == None:
        y_labels = ""
    seaborn.heatmap(xy_data, vmin=3.0, vmax=50.0, square=True, ax=axes, 
                    xticklabels=x_labels, yticklabels=y_labels, 
                    cmap="Blues", linewidth=0.1, linecolor="#f0f0ff", 
                    cbar_kws={"aspect": 40})
    cbar_axes = figure.axes[-1]
    cbar_axes.tick_params(labelsize=20)
    cbar_axes.set_xlabel("EF", size=24)
    pyplot.title("Enrichment Factor Heatmap", fontsize=24, y=1.01)
    pyplot.savefig(out_fn, bbox_inches="tight", pad_inches=0.3)


def read_ef(ef_csv):
    ef_dict = {}
    header = ef_csv.next()
    logging.info("Skipping input header: %s" % str(header))
    for row in ef_csv:
        tid, tname, event, ef, pvalue, qvalue = row
        ef = float(ef)
        if float(qvalue) > QVALUE_CUTOFF:
            continue
        if ef < EF_CUTOFF:
            continue
        if not tname.endswith("_HUMAN"):
            continue
        ef_dict[(event, tname)] = ef
    return ef_dict 

def normalize_ef(ef_dict):
    pass

def read_sea(sea_csv):
    pass

def normalize_sea(sea_dict):
    pass

def dual_heatmap(ef_csv, sea_csv, out_fn, alphabetical=False):
    """Plot heatmap of EF data."""
    value_dict = read_ef(ef_csv)
    event_ids = list(set(event for event, tname in value_dict.iterkeys()))
    logging.info("Read in %d event ids" % len(event_ids))
    target_ids = list(set(tname for event, tname in value_dict.iterkeys()))
    target_ids.sort()
    logging.info("Read in %d target ids" % len(target_ids))
    evalue_data = build_array(value_dict, event_ids, target_ids)
    clustered_event_ids = hierarchical_x(evalue_data, event_ids, 
                                         threshold=0.5)
    clustered_target_ids = hierarchical_y(evalue_data, target_ids, 
                                          threshold=0.9)
    graph_target_ids = clustered_target_ids
    if alphabetical:
        graph_target_ids = target_ids
    evalue_data = build_array(value_dict, clustered_event_ids,
                              graph_target_ids)
    plot_heatmap(evalue_data, out_fn, x_labels=clustered_event_ids,
                 y_labels=graph_target_ids)


def handler(ef_fn, sea_fn, out_fn, **kwargs):
    """I/O handling for the script."""
    ef_f = gopen(ef_fn)
    ef_csv = csv.reader(ef_f)
    sea_f = gopen(sea_fn)
    sea_csv = csv.reader(sea_f)
    try:
        try:
            dual_heatmap(ef_csv, sea_csv, out_fn, **kwargs)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        in_f.close()
    return 0


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Plot heatmap for EF analysis."
    parser = ArgumentParser(description=description)
    parser.add_argument("ef_file", 
                        help="input EF analysis CSV file")
    parser.add_argument("sea_file", 
                        help="input SEA results CSV file")
    parser.add_argument("outfile", 
                        help="output heatmap PNG file")
    parser.add_argument("-a", "--alphabetical", action="store_true",  
                        help="Arrange targets in alphabetical " +
                             "instead of clustured order")
    options = parser.parse_args(args=argv[1:])
    return handler(ef_fn=options.ef_file, sea_fn=options.sea_file, 
                   out_fn=options.outfile, 
                   alphabetical=options.alphabetical)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

