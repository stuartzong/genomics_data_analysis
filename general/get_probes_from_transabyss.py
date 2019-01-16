#! /usr/bin/env python

import os
import os.path
import datetime
import argparse
import csv
import logging

import colorlog
logger = colorlog.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter('%(log_color)s%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(handler)


def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter variants based on qulaity and somatic filters')
    parser.add_argument(
        '-i', '--input_file',
        help='specify input file', required=True)
    args = parser.parse_args()
    return args

def convert_probes(infile):
    outfile = '.'.join([infile, 'probe'])
    with open(outfile,  'wb') as writer:
        # writer = csv.writer(opf, delimiter='\t')
        with open(infile, 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    writer.write(line)
                else: # sequence
                    midpoint = len(line)/2
                    print line, midpoint, type(line)
                    print midpoint-25
                    probe = line[midpoint-25: midpoint+25]
                    writer.write(probe)
                    writer.write('\n')
 
def main():
    start = datetime.datetime.now()
    logger.info("script starts at: %s" % start)
    args = parse_args()
    transabyss_probe_file = args.input_file

    convert_probes(transabyss_probe_file)
    end = datetime.datetime.now()
    logger.info("script ends at: %s\n" % end)


if __name__ == '__main__':
    main()
