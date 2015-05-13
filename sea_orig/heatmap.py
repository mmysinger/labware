#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Generate a heatmap from original SEA scores file

Adapted from original code by Leo Gendelev

Michael Mysinger 201504 Created
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


module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen
from sea_orig.clustering import hierarchical_x, hierarchical_y, hierarchical_both

def build_array(value_dict, x_labels, y_labels):
    values = []
    for y_label in y_labels:
        values.append([value_dict.get((x_label, y_label), 0.0)
                            for x_label in x_labels])
    return numpy.array(values)

def plot_heatmap(xy_data, out_fn, x_labels=None, y_labels=None):
    logging.info("Input data shape: %s" % str(xy_data.shape))
    figure, axis = pyplot.subplots()
    figure.set_size_inches(10, 15)
    cmap = matplotlib.cm.get_cmap('Blues')
    cmap.set_under('w')
    heatmap = axis.imshow(xy_data, cmap=cmap, interpolation='None',
                          vmin=10.0, vmax=90.0, aspect=0.5)
    cbar = pyplot.colorbar(heatmap, aspect=40)
    cbar.set_label('$\mathregular{-Log_{10}}$ E-Value', size=24)
    cbar.ax.tick_params(labelsize=24) 
    axis.set_xticks(numpy.arange(xy_data.shape[1]), minor=False)
    axis.set_yticks(numpy.arange(xy_data.shape[0]), minor=False)
    if x_labels == None:
        x_labels = ""
    axis.set_xticklabels(x_labels, minor=False, rotation='vertical', fontsize=8)
    if y_labels == None:
        y_labels = ""
    axis.set_yticklabels(y_labels, minor=False, fontsize=8)
    pyplot.title('E-value Heatmap', fontsize=32)
    #axis.axis('tight')
    pyplot.savefig(out_fn)
    #pyplot.show()


def heatmap(in_csv, out_fn, ref_suffix="_10000"):
    """Module code starts here."""
    log_evalue_dict = {}
    header = in_csv.next()
    logging.info("Skipping input header: %s" % str(header))
    for row in in_csv:
        query_id, ref_id, evalue, maxtc = row
        if ref_suffix and ref_id.endswith(ref_suffix):
            ref_id = ref_id[:-len(ref_suffix)]
        log_evalue_dict[(query_id, ref_id)] = -numpy.log10(float(evalue))
    query_ids = list(set(query_id for query_id, ref_id in
                         log_evalue_dict.iterkeys()))
    logging.info("Read in %d query ids" % len(query_ids))
    ref_ids = list(set(ref_id for query_id, ref_id in
                       log_evalue_dict.iterkeys()))
    #ref_ids.sort()
    logging.info("Read in %d reference ids" % len(ref_ids))
    evalue_data = build_array(log_evalue_dict, query_ids, ref_ids)
    clustered_query_ids = hierarchical_x(evalue_data, query_ids, threshold=0.5)
    clustered_ref_ids = hierarchical_y(evalue_data, ref_ids, threshold=0.9)
    evalue_data = build_array(log_evalue_dict, clustered_query_ids,
                              clustered_ref_ids)
    plot_heatmap(evalue_data, out_fn, x_labels=clustered_query_ids,
                 y_labels=clustered_ref_ids)


def handler(in_fn=None, out_fn=None, **kwargs):
    """I/O handling for the script."""
    if in_fn is None:
        in_f = sys.stdin
    else:
        in_f = gopen(in_fn, "r")
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


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Template: paragraph description goes here."
    parser = ArgumentParser(description=description)
    parser.add_argument("-i", "--infile", default=None, 
                        help="input file (default: stdin)")
    parser.add_argument("-o", "--outfile", default=None, 
                        help="output file (default: stdout)")
    options = parser.parse_args(args=argv[1:])
    return handler(in_fn=options.infile, out_fn=options.outfile)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

