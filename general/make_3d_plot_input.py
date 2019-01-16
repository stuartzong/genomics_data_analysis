#! /usr/bin/env python

import os
import os.path
import datetime
import argparse
import csv
import logging
from pprint import pprint

import colorlog

import ConfigParser

logger = colorlog.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter('%(log_color)s%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(handler)


def write_files(variant_afs, outfile):
    with open(outfile, 'wb') as fh1:
        writer = csv.writer(fh1, delimiter='\t')
        headers = ['variant', 'diagnois_af', 'relapse_af', 'postmortem']
        writer.writerow(headers)
        for variant in variant_afs:
            if len(variant_afs[variant]) == 3:
                afs = [variant_afs[variant][patient_status]
                       for patient_status in variant_afs[variant]]
                print variant_afs[variant]
                print afs
                content = [variant] + [afs[-i] for i in range(1, 4)]
                print content
                writer.writerow(content)


def make_af_dict(infile):
    variant_afs = {}
    with open(infile) as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            gene = line['gene']
            patient_status = line['patientID']
            variant = '_'.join([line['gene'],
                                line["chromosome"],
                                line["position"],
                                line["referenceBase"],
                                line["alternativeBase"]])
            af = line['tumourDNAAlleleFrequency']
            try:
                variant_afs[variant][patient_status].append(af)
            except KeyError:
                if variant not in variant_afs:
                    variant_afs[variant] = {}
                if patient_status not in variant_afs[variant]:
                    variant_afs[variant][patient_status] = af
    pprint(variant_afs)
    return variant_afs


def parse_args():
    parser = argparse.ArgumentParser(
        description='Filter variants based on qulaity and somatic filters')
    parser.add_argument(
        '-i', '--input_file',
        help='specify input file which is unfiltered variant summary file',
        required=True)
    args = parser.parse_args()
    return args


def main():
    start = datetime.datetime.now()
    logger.info("Quality and somatic filtering script starts at: %s" % start)
    args = parse_args()
    filtered_summary = args.input_file
    outfile = '.'.join([filtered_summary, '3d'])
    variant_afs = make_af_dict(filtered_summary)
    write_files(variant_afs, outfile)

    
    end = datetime.datetime.now()
    logger.info("Quality and somatic filtering script ends at: %s\n" % end)


if __name__ == '__main__':
    main()
