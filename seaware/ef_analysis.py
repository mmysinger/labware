#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Compute enrichment factors and q-values.

Michael Keiser 201210 Created
Garrett Gaskins 201505 Converted p and q-value calculations to Python 
Michael Mysinger 201508 Refactored and speed optimized
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser

import csv
from collections import namedtuple, defaultdict, Counter
import numpy as np
from scipy import stats
import statsmodels.sandbox.stats.multicomp as sm_multi

module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen

Target = namedtuple("Target", "name description")

# XXX - MMM currently set for bitterdb, may want better defaults 
CUTOFF_EF = 3.0           # Nat2012: EF > 1
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
        cid, eid = row[:2]
        # XXX - may want two optional columns, for extra drug and event info
        # read Garrett's file format
        #cid, altid, eid = row
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
        # read Garrett's file format
        #cid, altid, tid, affinity, pvalue, maxtc, name = row
        #desc = ""
        if cid not in has_event:
            rejects.add(cid)
            continue
        # implicitly takes the union over remaining affinity groups
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
    """Prune event molecules that are not mapped to targets"""
    pruned = dict((event, drugs & has_target) for event, drugs in
                                                events_to_drugs.iteritems())
    rejects = has_event - flatten_setdict(pruned)
    logging.info("Pruned %d molecules that were not mapped to targets" % 
                 len(rejects))
    return pruned


def precompute_sums(events_to_drugs, targets_to_drugs):
    """Pre-compute E, T, and p sums for enrichment factor calculation"""
    # For each event, calculate total number of molecule-target pairs
    drugs_to_targets = flip_setdict(targets_to_drugs)
    drugs_to_events = flip_setdict(events_to_drugs)
    assert(len(drugs_to_events) == len(drugs_to_targets))
    E = {}
    for event, drugs in events_to_drugs.iteritems():
        E[event] = sum(len(drugs_to_targets[drug]) for drug in drugs)

    # For each target, calculate total number of molecule-event pairs
    T = {}
    for target, drugs in targets_to_drugs.iteritems():
        T[target] = sum(len(drugs_to_events[drug]) for drug in drugs)

    # For each target-event pair, calculate number of molecules in common
    p = {}
    for target, t_drugs in targets_to_drugs.iteritems():
        for event, e_drugs in events_to_drugs.iteritems():
            p[(target, event)] = len(t_drugs & e_drugs)    
    return E, T, p


def compute_efs(E, T, p, events_to_drugs, targets_to_drugs, 
                min_pairs=CUTOFF_MINPAIRS, ef_cutoff=CUTOFF_EF):
    """Compute enrichment factors"""
    # EF = p/(E*T/P)
    # P = global count of all event-molecule-target triplets
    # p = for each target-event pair, count common linking molecules 
    # E = for each event, number of linked molecule-target pairs
    # T = for each target, number of linked molecule-event pairs

    # Total count of all event-molecule-target triplets
    P = sum(p.values())
    logging.info("Total number event-molecule-target triplets: P = %d" % P)
    P = float(P)
    efs = {}
    for event in events_to_drugs:
        ee = E[event]
        if ee == 0:
            continue
        for target in targets_to_drugs:
            tt = T[target]
            pte = p[(target, event)] 
            if pte < min_pairs:
                continue
            try:
                ef = float(pte)/(ee*tt/P)
            except:
                logging.warn('EF Error', event, target, 
                             pte, ee, tt)
                continue
            efs[(target, event)] = ef
    logging.info("Computed %d target-event enrichment factors" % len(efs))
    return efs


def map_contingency_tables(efs, events_to_drugs, targets_to_drugs):
    """Calculate contingency table for every target-event pair"""
    # Count number of drug-target pairs for each drug
    target_counts = Counter()
    for drugs in targets_to_drugs.itervalues(): 
        target_counts.update(drugs)
    num_pairs = sum(target_counts.itervalues())
    # Use counts to quickly compute contingency sums
    logging.info("Computing contingency tables")
    contingency_tables = {}
    for target, event in efs.iterkeys():
        e_drugs = events_to_drugs[event]
        t_drugs = targets_to_drugs[target]
        both = sum(target_counts[x] for x in e_drugs & t_drugs) 
        events = sum(target_counts[x] for x in e_drugs) - both
        targets = sum(target_counts[x] for x in t_drugs) - both
        neither = num_pairs - both - events - targets
        contingency_tables[(target, event)] = np.array([[both, events], 
                                                        [targets, neither]])
    return contingency_tables


def compute_q_values(contingencies):
    """Compute p and q-values"""
    logging.info("Computing p and q-values")
    target_event_pairs = []
    p_vals = []
    for (target, event), table in contingencies.iteritems():
        chi2, pvalue, ddof, expected = stats.chi2_contingency(table)
        target_event_pairs.append((target, event))
        p_vals.append(pvalue)
    #Calculate the qvalue (p-adjusted FDR)
    reject_array, q_vals, alpha_c_sidak, alpha_c_bonf = sm_multi.multipletests(
        p_vals, alpha=0.05, method='holm')
    return target_event_pairs, p_vals, q_vals


def ef_analysis(events_reader, results_reader, min_pairs=CUTOFF_MINPAIRS, 
                ef_cutoff=CUTOFF_EF):
    """Compute enrichment factors and write q-values."""
    logging.info("Using min-pairs cutoff = %d" % min_pairs)
    logging.info("Using EF cutoff = %.2f" % ef_cutoff)
    events_to_drugs, has_event = read_events(events_reader)
    targets_to_drugs, has_target, targets = read_results(results_reader, 
                                                             has_event)
    events_to_drugs = prune_events(events_to_drugs, has_event, has_target)
    del has_target, has_event
    E, T, p = precompute_sums(events_to_drugs, targets_to_drugs)
    efs = compute_efs(E, T, p, events_to_drugs, targets_to_drugs, 
                      min_pairs=min_pairs, ef_cutoff=ef_cutoff)
    contingencies = map_contingency_tables(efs, events_to_drugs, 
                                           targets_to_drugs)
    target_event_pairs, p_vals, q_vals = compute_q_values(contingencies)
    assert(len(target_event_pairs) == len(efs))
    logging.info("Writing output")
    yield ["uniprot_id", "targ_name", "event", "ef", "p-value", "q-value"]
    count = 0
    for te_pair, p_val, q_val in zip(target_event_pairs, p_vals, q_vals):
        target, event = te_pair
        if efs[te_pair] > ef_cutoff:
            count += 1
            yield [target, targets[target].name, event, "%.5g" % efs[te_pair],
                   "%.5g" % p_val, "%.5g" % q_val]
    logging.info("Wrote %d q-values to output" % count)


def handler(events_fn, results_fn, out_fn=None, **kwargs):
    """I/O handling for the script."""
    logging.info("Events file: %s" % events_fn)
    events_f = open(events_fn, "r")
    events_reader = csv.reader(events_f)
    logging.info("SEAware results file: %s" % results_fn)
    results_f = open(results_fn, "r")
    results_reader = csv.reader(results_f)
    if out_fn is None:
        out_f = sys.stdout
        out_location = "standard out"
    else:
        out_f = open(out_fn, "w")
        out_location = out_fn
    logging.info("Output file: %s" % out_location)
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
    description = "Compute enrichment factors and q-values"
    parser = ArgumentParser(description=description)
    parser.add_argument("events",  
                        help="Events file mapping molecules to events")
    parser.add_argument("results",  
                        help="SEAware results mapping molecules to targets")
    parser.add_argument("-m", "--min-pairs", type=int, default=CUTOFF_MINPAIRS, 
        help="Minimum pairs cutoff for EF analysis (default: %(default)s)")
    parser.add_argument("-e", "--ef-cutoff", type=float, default=CUTOFF_EF, 
        help="EF factor cutoff to write results (default: %(default)s)")
    parser.add_argument("-o", "--outfile", default=None, 
                        help="output file (default: stdout)")
    options = parser.parse_args(args=argv[1:])
    return handler(events_fn=options.events, results_fn=options.results, 
                   out_fn=options.outfile, min_pairs=options.min_pairs, 
                   ef_cutoff=options.ef_cutoff)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

