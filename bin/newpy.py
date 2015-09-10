#!/usr/bin/env python
"""Generate an empty template of a python module.

Michael Mysinger 200707 Created
Michael Mysinger 201110 Make handleio more capable (exceptions, generator)
Michael Mysinger 201206 Modify for SeaChange
Michael Mysinger 201502 Modify for argparse
Michael Mysinger 201509 Switch to explicit file arguments, enhance logging
"""

import os
import sys
import time
import os.path as op
from argparse import ArgumentParser

OWNER = 'Michael Mysinger'

def create_template(filename, creator=OWNER, standalone=False):
    program_name = op.splitext(op.basename(filename))[0]
    if op.exists(filename):
        raise IOError("Error: file %s already exists!" % filename)
    today = time.strftime('%Y%m')

    shared_text = """
module_path = os.path.realpath(os.path.dirname(__file__)) 
labware_path = os.path.join(module_path, "..")
sys.path.append(labware_path)
from libraries.lab_utils import ScriptError, gopen
"""
    
    standalone_text = '''

class ScriptError(StandardError):
    def __init__(self, msg, value=99):
        self.value = value
        Exception.__init__(self, msg)


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
'''

    if standalone:
        util_text = standalone_text
    else:
        util_text = shared_text

    template = '''#!/usr/bin/env python
"""
Copyright (C) 2015 Michael M Mysinger

Template: one line summary goes here.

%s %s Created
"""

import os
import sys
import logging
import os.path as op
from argparse import ArgumentParser
import csv
%s

def %s(in_reader):
    """Module code starts here."""
    logging.info("Hello World!")
    for row in in_reader:
      yield row


def handler(in_fn=None, out_fn=None, **kwargs):
    """I/O handling for the script."""
    logging.info("Input file: %%s" %% in_fn)
    in_f = gopen(in_fn, "r")
    in_reader = csv.reader(in_f)
    logging.info("Output file: %%s" %% out_fn)
    out_f = gopen(out_fn, "w")
    out_writer = csv.writer(out_f)
    try:
        try:
            for result in %s(in_reader, **kwargs):
                out_writer.writerow(result)
        except ScriptError, message:
            logging.error(message)
            return message.value
    finally:
        in_f.close()
        out_f.close()
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
    log_format = "%%(levelname)s: %%(message)s"
    logging.basicConfig(level=log_level, format=log_format)
    description = "Template: paragraph description goes here."
    parser = ArgumentParser(description=description)
    parser.add_argument("infile", 
                        help="input CSV file")
    parser.add_argument("outfile", 
                        help="output CSV file")
    options = parser.parse_args(args=argv[1:])
    log_fn = options.outfile.replace(".csv", "") + ".log"
    add_file_logger(log_fn, level=log_level, format=log_format)
    return handler(in_fn=options.infile, out_fn=options.outfile)


if __name__ == "__main__":
    sys.exit(main(sys.argv))

''' % (creator, today, util_text, program_name, program_name)
    f = open(filename, 'w')
    f.write(template)
    f.close()
    # chmod +x filename
    mod = os.stat(filename)[0] & 0777 
    os.chmod(filename, mod | 0111)
    return True

def main(argv):
    """Parse arguments."""
    description = "Generate a python module template."
    parser = ArgumentParser(description=description)
    parser.add_argument("FILE", help="new Python script filename")
    parser.add_argument("-c", "--creator", default=OWNER, 
        help="set the name of the script creator (default: %(default)s)")
    parser.add_argument("-s", "--standalone", action="store_true", 
        help="Use standalone copy of the included utility functions")
    options = parser.parse_args(args=argv[1:])
    return not create_template(options.FILE, creator=options.creator,
                               standalone=options.standalone)

if __name__ == "__main__":
    sys.exit(main(sys.argv))
