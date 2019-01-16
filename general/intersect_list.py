#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice
import ConfigParser

def __main__():
    parser = argparse.ArgumentParser(description='Filter variants based on qulaity and somatic filters')
    parser.add_argument('-i1','--input_file1', help='specify input file', required=True)
    parser.add_argument('-i2','--input_file2', help='specify input file', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    infile_a = args['input_file1']
    infile_b = args['input_file2']

    a = get_list(infile_a)
    b = get_list(infile_b)

    shared = intersect(a, b)
    [uniq_a, uniq_b] = get_uniq(a, b)

    print "%s elements share are: \n%s" % (len(shared), shared)
    print "%s elements uniq to a is: \n%s" % (len(uniq_a), uniq_a)
    print "%s elements uniq to b is: \n%s" % (len(uniq_b), uniq_b)
 
def intersect(a, b):
    shared = list(set(a) & set(b))
    return shared

def get_uniq(a, b):
    uniq1= list(set(a) - set(b))
    uniq2 = list(set(b) - set(a))
    return [uniq1, uniq2]

def get_list(infile):
    a = []
    with open (infile, 'r') as fh:
         for line in fh:
             sl = line.split()
             a.append(sl[0])
    return a


if __name__ == '__main__':
    __main__()

