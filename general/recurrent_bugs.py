#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice


def parse_file(infile):
    bugs_dict = dict()
    with open (infile, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        header = records.fieldnames
        bugs = header[4:]
        for line in records:
            for bug in bugs:
                if (line[bug] != '0'):
                    patient = line['patient']
                    try:
                        bugs_dict[bug].append(patient)     
                    except KeyError:
                        bugs_dict[bug] = [patient]
        for bug in bugs:
            occurrence = len(list(set(bugs_dict[bug])))
            print "%s\t%s" % (occurrence, bug)




def __main__():
    print "BBT summary scripts starts at: %s\n" % datetime.datetime.now()     

    parser = argparse.ArgumentParser(description='generate BBT result matrix')
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    args = vars(parser.parse_args())

    # input file contains the full path to *.adjusted.normalized.txt files
    input_file = args['input_file']

    print "BBT input file is:\n%s\n" % (input_file)

    parse_file(input_file)

    print "Summarization scripts finished on: %s\n" % datetime.datetime.now()    

if __name__ == '__main__':
    __main__()

