#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice
import ConfigParser

def __main__():
    print "Quality and somatic filtering script starts at: %s\n" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='Filter variants based on qulaity and somatic filters')
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    parser.add_argument('-f','--somatic_file', help='specify possible somatic positions', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    variant_input_file = args['input_file']
    somatic_file = args['somatic_file']
    somatic_list = make_somatic_list(somatic_file)
    filter_variants(variant_input_file, somatic_list)


def make_somatic_list(somatic_file):
    somatic_list = []
    with open (somatic_file, 'r') as handle:
         for line in handle:
             sl = line.strip().split('\t')
             chr = sl[0]
             pos = sl[1]
             join = '_'.join([chr, pos])
             somatic_list.append(join)
    somatic_list = list(set(somatic_list))
    #print somatic_list
    return somatic_list

def filter_variants(filtered_summary, somatic_list):
    somatic_summary = ".".join([filtered_summary, "somatic" ])
    with open (somatic_summary,  'wb') as fh2:
        writer2 = csv.writer( fh2, delimiter='\t' )
        print "The filtered summary file is: %s.\n" % filtered_summary
        with open (filtered_summary, 'r') as handle:
             records = csv.DictReader(handle,  delimiter='\t')
             headers = records.fieldnames
             writer2.writerow(headers)
             for line in records:
                 gene = line['gene']
                 chr = line["chromosome"]
                 pos = line["position"]
                 ref = line["ref_base"]
                 alt = line["alt_base"]
                 var = '_'.join([chr, pos])
                 #print var
                 content = [line[i] for i in headers]
                 #quality filtering
                 if (var in somatic_list):
                     writer2.writerow(content)


if __name__ == '__main__':
    __main__()

