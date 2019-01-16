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
    parser.add_argument('-f','--patient_file', help='specify possible somatic positions', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    variant_input_file = args['input_file']
    patient_file = args['patient_file']
    patient_list = make_patient_list(patient_file)
    filter_variants(variant_input_file, patient_file,  patient_list)


def make_patient_list(patient_file):
    patient_list = []
    with open (patient_file, 'r') as handle:
         for line in handle:
             sl = line.strip().split('\t')
             pat = sl[0]
             status = sl[1]
             join = '_'.join([pat, status])
             patient_list.append(join)
    patient_list = list(set(patient_list))
    #print patient_list
    return patient_list

def filter_variants(filtered_summary, patient_file, patient_list):
    subset_summary = ".".join([filtered_summary, patient_file])
    with open (subset_summary,  'wb') as fh2:
        writer2 = csv.writer( fh2, delimiter='\t' )
        print "The filtered summary file is: %s.\n" % filtered_summary
        with open (filtered_summary, 'r') as handle:
             records = csv.DictReader(handle,  delimiter='\t')
             headers = records.fieldnames
             writer2.writerow(headers)
             for line in records:
                 pat = line['patient_ID']
                 #print var
                 content = [line[i] for i in headers]
                 #quality filtering
                 if (pat in patient_list):
                     writer2.writerow(content)


if __name__ == '__main__':
    __main__()

