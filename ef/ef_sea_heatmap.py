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
from sea_orig.heatmap import read_sea
from ef_heatmap import build_array, read_ef


def summed_array(value_dicts, x_labels, y_labels):
    values = []
    for y_label in y_labels:
        values.append([sum(value_dict.get((x_label, y_label), 0.0) 
            for value_dict in value_dicts) for x_label in x_labels])
    return numpy.array(values)


def joined_array(value_dicts, x_labels, y_labels):
    values = []
    for y_label in y_labels:
        y_values = []
        for x_label in x_labels:
            for value_dict in value_dicts:
                y_values.append(value_dict.get((x_label, y_label), 0.0))
        values.append(y_values)
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
    seaborn.heatmap(xy_data, vmin=-1.0, vmax=1.0, square=True, ax=axes, 
                    xticklabels=x_labels, yticklabels=y_labels, 
                    cmap="RdBu", linewidth=0.1, linecolor="#f0f0ff", 
                    cbar=False)
    pyplot.title("Dual SEA & EF Heatmap", fontsize=24, y=1.01)
    pyplot.savefig(out_fn, bbox_inches="tight", pad_inches=0.3)


def normalize(in_dict, lower=0.0, percentile=90.0):
    """Normalize dictionary values from 0 to 1 (-1 if inverted)"""
    upper = numpy.percentile(in_dict.values(), percentile)
    height = upper - lower
    out_norm = {}
    for k, v in in_dict.iteritems():
        if v >= upper:
            norm = 1.0
        elif v <= lower:
            norm = 0.0
        else:
            norm = (v - lower)/height
        out_norm[k] = norm
    return out_norm, (lower, upper)


def invert(in_dict):
    return dict((k, -v) for k, v in in_dict.iteritems())        


def dual_heatmap(ef_csv, sea_csv, out_fn, alphabetical=False, clusters="both"):
    """Plot heatmap of EF data.

    clusters = both | sea | ef
    """
    ef_dict = read_ef(ef_csv, clusters=False)
    ef_norm, limits = normalize(ef_dict, lower=3.0)
    logging.info("EFs normalized from %g to %g" % limits) 
    ef_events = set(event for event, target in ef_norm.iterkeys())
    ef_targets = set(target for event, target in ef_norm.iterkeys())
    logging.info("Found %d event and %d target labels in EF data" % 
                 (len(ef_events), len(ef_targets)))
    sea_dict = read_sea(sea_csv)
    sea_norm, limits = normalize(sea_dict, lower=6.0)
    logging.info("SEA -log10 p-values normalized from %g to %g" % limits) 
    sea_events = set(event for event, target in sea_norm.iterkeys())
    sea_targets = set(target for event, target in sea_norm.iterkeys())
    logging.info("Found %d event and %d target labels in SEA data" % 
                 (len(sea_events), len(sea_targets)))
    event_ids = list(ef_events | sea_events)
    target_ids = list(ef_targets | sea_targets)
    logging.info("Totaled %d event and %d target labels" % 
                 (len(event_ids), len(target_ids)))
    if clusters.lower() == "both":
        logging.info("Clustering using both EF and SEA results")
        value_data = summed_array([sea_norm, ef_norm], event_ids, target_ids)
    elif clusters.lower() == "sea":
        logging.info("Clustering using SEA results only")
        value_data = build_array(sea_dict, event_ids, target_ids)
    elif clusters.lower() == "ef":
        logging.info("Clustering using EF results only")
        value_data = build_array(ef_dict, event_ids, target_ids)
    else:
        raise ValueError("Unrecognized clusters parameter %s" % clusters)
    clustered_event_ids = hierarchical_x(value_data, event_ids, 
                                         threshold=0.5)
    clustered_target_ids = hierarchical_y(value_data, target_ids, 
                                          threshold=0.9)
    graph_target_ids = clustered_target_ids
    if alphabetical:
        logging.info("Showing alphabetical target list instead of clusters")
        target_ids.sort()
        graph_target_ids = target_ids
    ef_norm = invert(ef_norm)
    value_data = joined_array([sea_norm, ef_norm], 
                               clustered_event_ids, graph_target_ids)
    graph_event_ids = []
    for event in clustered_event_ids:
        graph_event_ids.append(event + " sea")
        graph_event_ids.append(event + " ef")
    plot_heatmap(value_data, out_fn, x_labels=graph_event_ids,
                 y_labels=graph_target_ids)


def handler(ef_fn, sea_fn, out_fn, **kwargs):
    """I/O handling for the script."""
    logging.info("EF input file: %s" % ef_fn)
    ef_f = gopen(ef_fn)
    ef_csv = csv.reader(ef_f)
    logging.info("SEA input file: %s" % sea_fn)
    sea_f = gopen(sea_fn)
    sea_csv = csv.reader(sea_f)
    logging.info("Output file: %s" % out_fn)
    try:
        try:
            dual_heatmap(ef_csv, sea_csv, out_fn, **kwargs)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        ef_f.close()
        sea_f.close()
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
    parser.add_argument("ef_file", 
                        help="EF analysis qvalues CSV file")
    parser.add_argument("sea_file", 
                        help="SEA scores CSV file")
    parser.add_argument("outfile", 
                        help="output heatmap PNG file")
    parser.add_argument("-a", "--alphabetical", action="store_true",  
                        help="Arrange targets in alphabetical " +
                             "instead of clustured order")
    options = parser.parse_args(args=argv[1:])
    log_fn = options.outfile.replace(".png", "") + ".log"
    add_file_logger(log_fn, level=log_level, format=log_format)
    return handler(ef_fn=options.ef_file, sea_fn=options.sea_file, 
                   out_fn=options.outfile, 
                   alphabetical=options.alphabetical)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

