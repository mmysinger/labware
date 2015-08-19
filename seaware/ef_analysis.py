#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Template: one line summary goes here.

Michael Mysinger 201505 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser
import csv
from collections import namedtuple, defaultdict

module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen

Target = namedtuple("Target", "name description")

# XXX - MMM currently set for bitterdb, may want better defaults 
CUTOFF_EF = 2.0           # Nat2012: EF > 1
CUTOFF_MINPAIRS = 3       # Nat2012: Target-ADR pairs > 10 retained


def flip_setdict(in_dict):
    """Reverse a dict of sets"""
    out_setdict = defaultdict(set)
    for k, v in in_dict.iteritems():
        for x in v:
            out_setdict[x].add(k)
    return out_setdict


def flatten_setdict(in_dict):
    """Flatten a dict of sets"""
    out_set = set()
    for x in in_dict.values():
        out_set.update(x)
    return out_set


def read_events(events_reader):
    """Read events to molecules mapping"""
    logging.info("Reading events")
    events_to_drugs = defaultdict(set)
    for row in events_reader:
        # XXX - Is affinity column needed?
        cid, eid, affinity, event = row
        events_to_drugs[eid].add(cid)
    has_event = flatten_setdict(events_to_drugs)
    logging.info("Mapped %d events to %d molecules" % (
        len(events_to_drugs), len(has_event)))
    return events_to_drugs, has_event


def read_results(results_reader, has_event):
    """Read targets to molecules mapping"""
    logging.info("Reading targets")
    header = results_reader.next()
    logging.info("Skipping SEAware results header: %s" % str(header))
    targets = {}
    targets_to_drugs = defaultdict(set)
    rejects = set()
    for row in results_reader:
        cid, smiles, tid, affinity, pvalue, maxtc, name, desc = row
        if cid not in has_event:
            rejects.add(cid)
            continue
        targets_to_drugs[tid].add(cid)
        if tid not in targets:
            targets[tid] = Target(name, desc)
    has_target = flatten_setdict(targets_to_drugs)
    logging.info("Mapped %d targets to %d molecules" % (
        len(targets_to_drugs), len(has_target)))
    logging.info("Skipped %d molecules that were not mapped to events" % 
                 len(rejects))
    return targets_to_drugs, has_target, targets


def prune_events(events_to_drugs, has_event, has_target):
    pruned = dict((event, drugs & has_target) for event, drugs in
                                                events_to_drugs.iteritems())
    rejects = has_event - flatten_setdict(pruned)
    logging.info("Pruned %d molecules that were not mapped to targets" % 
                 len(rejects))
    return pruned


def precompute_stats(events_to_drugs, targets_to_drugs):
    # For each event, calculate total number of molecule-target pairs
    drugs_to_targets = flip_setdict(targets_to_drugs)
    E = {}
    for event, drugs in events_to_drugs.iteritems():
        E[event] = sum(len(drugs_to_targets[drug]) for drug in drugs)

    # For each target, calculate total number of molecule-event pairs
    drugs_to_events = flip_setdict(events_to_drugs)
    assert(len(drugs_to_events) == len(drugs_to_targets))
    del drugs_to_targets
    T = {}
    for target, drugs in targets_to_drugs.iteritems():
        T[target] = sum(len(drugs_to_events[drug]) for drug in drugs)
    del drugs_to_events

    # For each target-event pair, calculate number of molecules in common
    p = {}
    for target, t_drugs in targets_to_drugs.iteritems():
        for event, e_drugs in events_to_drugs.iteritems():
            p[(target, event)] = len(t_drugs & e_drugs)    
    return E, T, p

def ef_analysis(events_reader, results_reader, min_pairs=CUTOFF_MINPAIRS, 
                ef=CUTOFF_EF):
    """Enrichment factor analysis."""
    events_to_drugs, has_event = read_events(events_reader)
    targets_to_drugs, has_target, targets = read_results(results_reader, 
                                                             has_event)
    events_to_drugs = prune_events(events_to_drugs, has_event, has_target)
    del has_target, has_event
    E, T, p = precompute_stats(events_to_drugs, targets_to_drugs)
    P = sum(p.values())
    logging.info("Total number of molecular links between " + 
                 "all target-event pairs: P = %d" % P)
    pass


def handler(events_fn, results_fn, out_fn=None, **kwargs):
    """I/O handling for the script."""
    logging.info("Events file: %s" % events_fn)
    events_f = open(events_fn, "r")
    events_reader = csv.reader(events_f)
    logging.info("SEAware results file: %s" % results_fn)
    results_f = open(results_fn, "r")
    results_reader = csv.reader(results_f)
    if outfile is None:
        out_f = sys.stdout
        out_location = "standard out"
    else:
        out_f = open(outfile, "w")
        out_location = out_fn
    logging.info("Output file: %s" %% out_location)
    out_writer = csv.writer(out_f)
    try:
        try:
            for result in ef_analysis(events_reader, results_reader, 
                                      **kwargs):
                out_writer.writerow(result)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        events_f.close()
        results_f.close()
        out_f.close()
    return 0


def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Compute enrichment factors."
    parser = ArgumentParser(description=description)
    parser.add_argument("events",  
                        help="Events file mapping molecules to events")
    parser.add_argument("results",  
                        help="SEAware results mapping molecules to targets")
    parser.add_argument("-m", "--min-pairs", type=int, default=CUTOFF_MINPAIRS, 
                        help="Minimum pairs cutoff for EF analysis.")
    parser.add_argument("-e", "--ef", type=float, default=CUTOFF_EF, 
                        help="EF factor cutoff to write results.")
    parser.add_argument("-o", "--outfile", default=None, 
                        help="output file (default: stdout)")
    options = parser.parse_args(args=argv[1:])
    return handler(events_fn=options.events, results_fn=options.results, 
                   out_fn=options.outfile, min_pairs=options.min_pairs, 
                   ef=options.ef)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

