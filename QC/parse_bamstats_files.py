#!/usr/bin/env python

import datetime
import csv
import argparse
import re
import sys


def parse_bamstats(bamstats_files, bamstats_summary):
    with open(bamstats_summary,  'w') as opf:
        writer = csv.writer(opf,  delimiter='\t')
        headers = ['lib', 'aligned_bases', 'coverage', 'total_bases',
                   'Read_length', 'alignment_rate', 'duplication_rate']
        # print(type(headers))
        writer.writerow(headers)
        headers = '\t'.join(headers)
        print(headers)
        with open(bamstats_files, 'r') as fh:
            for f in fh:
                bamstats = f.strip()
                sl = bamstats.split('/')
                lib = sl[6]
                if re.match('^A[0-9][0-9][0-9][0-9][0-9]$', lib):
                    print("matched {0}".format(lib))
                else:
                    print("not matched {0}".format(lib))
                    lib = [i for i in sl if re.match('^A[0-9][0-9][0-9][0-9][0-9]', i)][0]
                    # lib = [i for i in sl if re.match('^A[0-9][0-9][0-9][0-9][0-9]$', i)][0]
                    # sys.exit()
                    print('rematched lib {0}'.format(lib))
                # lib = bamstats.split("/")[4]
                with open(bamstats,  'r') as handle:
                    for line in handle:
                        # print line
                        if line.startswith("Read_length:"):
                            read_length = int(line.strip().
                                              split("\t")[1].split(":")[0])
                        elif line.startswith("Total_Number_Of_Reads:"):
                            total_reads = int(line.strip().split("\t")[1])
                            # print "total_reads: %s" % total_reads
                        elif line.startswith("Number_Reads_Aligned:"):
                            aligned_reads = int(line.strip().split("\t")[1])
                            # print "aligned_reads: %s" % aligned_reads
                        elif line.startswith("Number_of_Duplicates:"):
                            dup_reads = int(line.strip().split("\t")[1])
                        # elif line.startswith("Estimate_for_X_coverage:"):
                        elif line.startswith("Estimate_for_genome_X_coverage:"):
                            coverage = "{:.1f}".format(float(line.strip().
                                                             split("\t")[1]))
                        elif line.startswith("Estimate_for_X_coverage:"):
                            coverage  = "{:.1f}".format(float(line.strip().
                                                              split("\t")[1]))
                        else:
                            continue
                    aligned_bases = "{:.5f}".format(aligned_reads
                                                    * read_length/1000000000.0)
                    total_bases = "{:.5f}".format(total_reads
                                                  * read_length/1000000000.0)
                    dup_rate = "{:.3f}".format(1.00* dup_reads/total_reads)
                    alignment_rate = "{:.3f}".format(1.00
                                                     * aligned_reads/total_reads)
                    content = [lib, aligned_bases, coverage, total_bases,
                               str(read_length), alignment_rate, dup_rate]
                    print('\t'.join(content))
                    writer.writerow(content)
      
def parse_args():
    parser = argparse.ArgumentParser(
        description='Parse DNA bamstats to report coverage, duplicate rate etc')
    parser.add_argument(
        '-i', '--input_file',
        help='specify input file, which contains bamstats files',
        required=True)
    parser.add_argument(
        '-m', '--message',
        help='This reports non-duplicate coverage extracted from bamstats file.',
        required=False)
    args = parser.parse_args()
    return args


def __main__():
    print("Scripts starts at: {0}!".format(datetime.datetime.now()))
    args = parse_args()
    bamstats_files = args.input_file
    bamstats_summary = "bamstats_summary.txt"
    parse_bamstats(bamstats_files, bamstats_summary)
    print("Scripts ends at: {0}!".format(datetime.datetime.now()))


if __name__ == '__main__':
    __main__()

