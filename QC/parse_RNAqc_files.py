#!/usr/bin/env python


import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from pprint import pprint
from itertools import islice
import operator
#import ConfigParser
from collections import defaultdict


def get_rnaqc_metrics(rnaqc_files, rnaqc_summary):
    with open(rnaqc_summary,  'w') as opf:
        writer = csv.writer(opf,  delimiter='\t')
        new_headers = ["library_name",
                       "Aligned_Gbp",
                       "total_gbp",
                       "exon_intron_ratio",
                       "num_genes_coverage_ge_10X"]
        writer.writerow(new_headers)
        joined_header = "\t".join(new_headers)
        with open(rnaqc_files, 'r') as fh:
            for line in fh:
                rnaqc_file = line.strip()
                library = rnaqc_file.split("/")[4]
                (library_name,
                 exon_intron_ratio,
                 num_gene_10x,
                 chastity_passed_gbps,
                 percent_chastity_passed,
                 total_gbps,
                 percent_aligned,
                 aligned_gbps,
                 headers) = parse_rnaqc_file(rnaqc_file)
                if (library == library_name):
                    print("\t".join([library,
                                     str(aligned_gbps),
                                     str(total_gbps),
                                     exon_intron_ratio,
                                     num_gene_10x]))
                    writer.writerow([library,
                                     aligned_gbps,
                                     total_gbps,
                                     exon_intron_ratio,
                                     num_gene_10x])
    print("available qc matrix include:")
    for i in headers:
        print("\t".join([str(headers.index(i) + 1), i]))


def parse_rnaqc_file(rnaqc_file):
    with open(rnaqc_file,  'r') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        headers = records.fieldnames
        for line in records:
            library_name = line[headers[0]]
            exon_intron_ratio = line[headers[19]]
            num_gene_10x = line[headers[16]]
            chastity_passed_gbps = float(line[headers[2]])
            percent_chastity_passed = float(line[headers[3]])
            total_gbps = "{:.2f}".format(100 *
                                         chastity_passed_gbps /
                                         percent_chastity_passed)
            percent_aligned = float(line[headers[4]])
            aligned_gbps = "{:.2f}".format(chastity_passed_gbps *
                                           percent_aligned/100)
    return [library_name,
            exon_intron_ratio,
            num_gene_10x,
            chastity_passed_gbps,
            percent_chastity_passed,
            total_gbps,
            percent_aligned,
            aligned_gbps,
            headers]


def parse_args():
    parser = argparse.ArgumentParser(
        description='Parse rnaqc file to report qc metrics!')
    parser.add_argument(
        '-i', '--input_file',
        help='specify input file which list all rnaqc file paths',
        required=True)
    args = parser.parse_args()
    return args


def __main__():
    print("Scripts starts at: {0}!".format(datetime.datetime.now()))
    args = parse_args()
    rnaqc_files = args.input_file
    rnaqc_summary = "rnaqc_summary.txt"
    get_rnaqc_metrics(rnaqc_files, rnaqc_summary)
    print("Scripts ends at: {0}!".format(datetime.datetime.now()))


if __name__ == '__main__':
    __main__()

