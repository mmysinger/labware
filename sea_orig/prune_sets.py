#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Prune SEA original inputs to common subset of molecule ids 
    (between smi, fp, and set files)

Michael Mysinger 201504 Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

from sea.util.load_util import buildFilename
from sea.util.util import *

class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)

# New gopen function
def gopen(fname, mode='rU'):
    """Flexibly open compressed files."""
    if os.path.splitext(fname)[-1] == '.%s' % 'gz':
        if 'b' not in mode:
            # gzip does not support universal newlines
            mode = mode.replace('U', '') + 'b'
        f = gzip.open(fname, mode)
    elif os.path.splitext(fname)[-1] == '.%s' % 'bz2':
        if 'b' not in mode:
            # bzip2 does support universal newlines
            mode = mode + 'b'
        f = bz2.BZ2File(fname, mode)
    else:
        f = open(fname, mode)
    return f

def prune_sets(set_filename, minimum_size=1):
    """Prune SEA original inputs to common subset of molecule ids"""
    # two pass
    # first pass: read ids

    in_set_fn = buildFilename(set_filename, SUFFIX_SET)
    out_set_fn =  in_set_fn.replace('.%s' % SUFFIX_SET,
                                    '_pruned.%s' % SUFFIX_SET)
    
    smiles_ids = set()
    in_fn = buildFilename(in_set_fn, SUFFIX_SMI)
    with gopen(in_fn) as in_f:
        for line in in_f:
            smi, cid = line.strip().split(SEP)
            smiles_ids.add(cid)

    fingerprint_ids = set()
    in_fn = buildFilename(in_set_fn, SUFFIX_FP)
    with gopen(in_fn) as in_f:
        for line in in_f:
            fp, cid = line.strip().split(SEP)
            fingerprint_ids.add(cid)            

    set_ids = set()
    in_fn = buildFilename(in_set_fn, SUFFIX_SET)
    with gopen(in_fn) as in_f:
        for line in in_f:
            tid, name, cids = line.strip().split(SEP)
            set_ids.update(cids.split(SUBSEP))

    # second pass: filter ids
    common_ids = smiles_ids & fingerprint_ids & set_ids

    in_fn = buildFilename(in_set_fn, SUFFIX_SMI)
    out_fn = out_set_fn.replace(SUFFIX_SET, SUFFIX_SMI)
    with gopen(in_fn) as in_f, gopen(out_fn, 'w') as out_f:
        for line in in_f:
            smi, cid = line.strip().split(SEP)
            if cid in common_ids:
                out_f.write(line)

    in_fn = buildFilename(in_set_fn, SUFFIX_FP)
    out_fn = out_set_fn.replace(SUFFIX_SET, SUFFIX_FP)
    with gopen(in_fn) as in_f, gopen(out_fn, 'w') as out_f:
        for line in in_f:
            fp, cid = line.strip().split(SEP)
            if cid in common_ids:
                out_f.write(line)

    in_fn = in_set_fn
    out_fn = out_set_fn
    with gopen(in_fn) as in_f, gopen(out_fn, 'w') as out_f:
        for line in in_f:
            tid, name, cids = line.strip().split(SEP)
            filtered_cids = [cid for cid in cids.split(SUBSEP)
                             if cid in common_ids]
            if len(filtered_cids) >= minimum_size:
                out_cids = SUBSEP.join(filtered_cids)
                out_f.write('%s\n' % SEP.join([tid, name, out_cids]))

def main(argv):
    """Parse arguments."""
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s: %(message)s")
    description = "Prune SEA original inputs to common subset of molecule ids"
    parser = ArgumentParser(description=description,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("set_filename", help="SEA .set filename")
    parser.add_argument("-m", "--minimum-size", type=int, default=1, 
                        help="minimum set size for new set file")
    options = parser.parse_args(args=argv[1:])
    return prune_sets(options.set_filename, minimum_size=options.minimum_size)

if __name__ == "__main__":
    sys.exit(main(sys.argv))

