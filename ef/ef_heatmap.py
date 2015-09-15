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

# Add labware to sys.path
module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)

from libraries.lab_utils import ScriptError, gopen
from sea_orig.clustering import hierarchical_x, hierarchical_y, \
    hierarchical_both
from sea_orig.heatmap import build_array

QVALUE_CUTOFF = 1.0e-3
EF_CUTOFF = 5.0


def build_array(value_dict, x_labels, y_labels):
    values = []
    for y_label in y_labels:
        values.append([value_dict.get((x_label, y_label), 0.0)
                            for x_label in x_labels])
    return numpy.array(values)


def read_ef(ef_csv, ef_cutoff=EF_CUTOFF, qvalue_cutoff=QVALUE_CUTOFF, 
            clusters=True, human=True):
    ef_dict = {}
    header = ef_csv.next()
    logging.info("Skipping EF input header: %s" % str(header))
    for row in ef_csv:
        tid, tname, event, ef, pvalue, qvalue = row
        ef = float(ef)
        if float(qvalue) > qvalue_cutoff:
            continue
        if ef < ef_cutoff:
            continue
        if human and not tname.endswith("_HUMAN"):
            continue
        if not clusters and event.startswith("cluster_"):
            continue
        ef_dict[(event, tname)] = ef
    logging.info("Read %d enrichment factors" % len(ef_dict))
    return ef_dict 


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
    seaborn.heatmap(xy_data, vmin=3.0, vmax=40.0, square=True, ax=axes, 
                    xticklabels=x_labels, yticklabels=y_labels, 
                    cmap="Blues", linewidth=0.1, linecolor="#f0f0ff", 
                    cbar_kws={"aspect": 40})
    cbar_axes = figure.axes[-1]
    cbar_axes.tick_params(labelsize=20)
    cbar_axes.set_xlabel("EF", size=24)
    pyplot.title("Enrichment Factor Heatmap", fontsize=24, y=1.01)
    pyplot.savefig(out_fn, bbox_inches="tight", pad_inches=0.3)


def heatmap(in_csv, out_fn, alphabetical=False):
    """Plot heatmap of EF data."""
    value_dict = read_ef(in_csv)
    event_ids = list(set(event for event, tname in value_dict.iterkeys()))
    target_ids = list(set(tname for event, tname in value_dict.iterkeys()))
    logging.info("Found %d event and %d target labels in EF data" % 
                 (len(event_ids), len(target_ids)))
    value_data = build_array(value_dict, event_ids, target_ids)
    clustered_event_ids = hierarchical_x(value_data, event_ids, 
                                         threshold=0.5)
    clustered_target_ids = hierarchical_y(value_data, target_ids, 
                                          threshold=0.9)
    graph_target_ids = clustered_target_ids
    if alphabetical:
        logging.info("Showing alphabetical target list instead of clusters")
        target_ids.sort()
        graph_target_ids = target_ids
    value_data = build_array(value_dict, clustered_event_ids,
                              graph_target_ids)
    plot_heatmap(value_data, out_fn, x_labels=clustered_event_ids,
                 y_labels=graph_target_ids)


def handler(in_fn=None, out_fn=None, **kwargs):
    """I/O handling for the script."""
    in_f = gopen(in_fn)
    in_csv = csv.reader(in_f)
    try:
        try:
            heatmap(in_csv, out_fn, **kwargs)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        in_f.close()
    return 0


def add_file_logger(filename, level=logging.INFO,
                    format="%%(levelname)s: %%(message)s",
                    logger=logging.getLogger()):
    file_handler = logging.FileHandler(filename, mode="w")
    log_formatter = logging.Formatter(format)
    file_handler.setFormatter(log_formatter)
    file_handler.setLevel(level)
    logger.addHandler(file_handler)


def main(argv):
    """Parse arguments."""
    log_level = logging.INFO
    log_format = "%(levelname)s: %(message)s"
    logging.basicConfig(level=log_level, format=log_format)
    description = "Plot heatmap for EF analysis."
    parser = ArgumentParser(description=description)
    parser.add_argument("infile", 
                        help="input EF analysis CSV file")
    parser.add_argument("outfile", 
                        help="output heatmap PNG file")
    parser.add_argument("-a", "--alphabetical", action="store_true",  
                        help="Arrange targets in alphabetical " +
                             "instead of clustured order")
    options = parser.parse_args(args=argv[1:])
    log_fn = options.outfile.replace(".png", "") + ".log"
    add_file_logger(log_fn, level=log_level, format=log_format)
    return handler(in_fn=options.infile, out_fn=options.outfile, 
                   alphabetical=options.alphabetical)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

