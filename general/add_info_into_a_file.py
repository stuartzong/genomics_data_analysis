#!/usr/bin/env python

import os
import stat
import os.path
import time
import datetime
import subprocess
import re
import sys
import glob
import argparse
import csv
from collections import defaultdict
from pprint import pprint
from itertools import islice
import fileinput
import shutil

# dict: key -> info to be added 
def make_key_info_dict(infile):
    print "infile is: %s" % infile
    with open(infile) as handle:
        info_dict = dict()
        records = csv.DictReader(handle,  delimiter='\t')
        headers = records.fieldnames
        # assume column 1 and 2 are the key and info
        print "xxxx", headers
        for line in records:
                #status = line[headers[4]]
                #if (status.lower() != 'normal'):
                #sl = line.split('\t')
                # key = '_'.join([line[headers[0]], line[headers[1]]]).replace(' ', '')
                # print "yyyyyyyyyy"
                key= line[headers[0]]
                  
                #info = [line[headers[0]],  line[headers[1]]]
                #info = '_'.join([line[headers[1]], line[headers[2]]])
                # info = "\t".join([line[headers[0]], line[headers[1]]])
                # info = line['Human_papillomavirus']
                info = line[headers[1]]
                info_dict[key] = info
                #info_dict[key2] = info
    print "There are %s key and value pairs in input file: %s." % (len(info_dict), infile)
    print "The key-info dict is: "
    pprint(info_dict)
    return info_dict

def add_info_to_file(infile, info_dict):
    outfile = ".".join([infile, "added"])
    writer = csv.writer(open(outfile, 'wb'),delimiter='\t')
    with open(infile) as handle:
        # records = csv.DictReader(filter(lambda row: row[0]!='#', handle), delimiter='\t')
        records = csv.DictReader(handle,  delimiter='\t')
        headers = records.fieldnames
        # new_headers = headers[:6] + headers[7:] + ['gene']
        new_headers = headers + ["added_info"]
        writer.writerow(new_headers)
        # assume id is the first column
        for line in records:
            #sl = line.split('\t')
            #key = line[headers[11]].split("_")[0]
            # key = '_'.join([line[headers[0]], line[headers[1]]]).replace(' ', '')
            # key = line[headers[0]].replace(' ', '')
            # keys = line[headers[2]].split(';')
            #keys = [key.replace(' ', '') for key in keys]
            #info = [";".join([info_dict[key] for key in keys])]
            # keys = line[headers[4]].split(',')
            # keys = [i.split('/')[0] for i in keys if i.startswith('NM_')]
            # key = line[header[4]]
            # try:
            #     info = ','.join(list(set([info_dict[i] for i in keys]))) 
            # except KeyError:
            #     info = ','.join(keys)
            # key1 = line['DNA_lib']
            # key2 = line['RNA_lib']
   
            # info1 = info_dict[key1]
            # info2 = info_dict[key2]
            # info = [info1] + [info2]
            key = line[headers[4]]
            info =info_dict[key]






            """
            type = line[headers[3]]
            if (type == 'transcriptome'):
                info = info_dict[key][2]
            elif (type == 'genome'):
                info = info_dict[key][1]
            """
            """ 
            print "aa", key
            try:
                info = info_dict[key]
            except KeyError:
                info = "FLT3-ITD_data_unavailable"
            print "bb", info
            """
            # ATTENTION: headers must be all unique, if two columns have the same name, this outputs wrong content 
            content = [line[i] for i in headers] + [info] 
            # content = [line[i] for i in headers] + [info1] + [info2]
            writer.writerow(content)

def swap_names(infile,dict_name):
  with open (infile) as reader:
    for line in reader:
        # Split the line along whitespace
        # Note: this fails if your filenames have whitespace
       line=line.rstrip()
       line=line.split()
       dict_key=line[0]
       line_list= line[1].split(',')
       print line_list
       dict_name[dict_key] = line_list
    print dict_name


#mapping = {}
#populate_dict('old_new_name.tsv',mapping)
#be_renamed=dict()
##populate_dict('file2modify.tsv',be_renamed)
#
#swap_names('file2modify.tsv',be_renamed)
#new_dict=dict()
#for key in be_renamed.keys():
#  aa_list=[]
#  for element in be_renamed[key]:
#    if element not in mapping.keys():
#      aa_list.append(element)
#      #print key,"\t",element 
#    if element in mapping.keys():
#      aa_list.append(mapping[element]) 
#     
#      #print key,"\t",mapping[element]
#  new_dict[key]=aa_list
#  #print key,"\t",aa_list
#for key in sorted(new_dict.keys()):
#  print key,"\t",new_dict[key]

'''# List the files in the current directory
for filename in os.listdir('.'):
    print filename
    root, extension = os.path.splitext(filename)
    if not root.startswith(prefix):
        # File doesn't end with this suffix; ignore it
        continue
    # Strip off the number of characters that make up suffix
    stripped_root = root[len(prefix):]
    if stripped_root in mapping:
        os.rename(filename, ''.join(mapping[stripped_root] + extension))

'''

def add_index_to_file(infile):
    outfile = ".".join([infile, "added"])
    writer = csv.writer(open(outfile, 'wb'),delimiter='\t')
    genes = dict()
    with open(infile) as handle:
        for line in handle:
            sl = line.split('\t')
            print sl
            key = sl[3].strip()
            info = "\t".join(sl)
            try:
                genes[key].append(info) 
            except KeyError:
                genes[key] = [info]
                
            # ATTENTION: headers must be all unique, if two columns have the same name, this outputs wrong content 
            # content = [line[i] for i in headers] + [info] 
            # content = [line[i] for i in headers] + [info1] + [info2]
            # writer.writerow(content)
        for gene in genes:
            intervals = genes[gene]
            for interval in intervals:
                index = intervals.index(interval)+1
                interval = interval.split('\t')
                content = interval[:3] + ['_'.join([gene, "EXON", str(index)])]
                writer.writerow(content)
        


def __main__():
    print "Scripts starts at: %s!" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='Add more info into a file based on a line identifier, eg. sample id')
    parser.add_argument('-i1','--input_file1', help='specify input file1', required=True)
    parser.add_argument('-i2','--input_file2', help='specify input file2', required=True)
    args = vars(parser.parse_args())

    input_file1 = args['input_file1']
    input_file2 = args['input_file2']

        
    #patient_clinics_dict = dict()
    patient_clinics_dict = make_key_info_dict(input_file1)

    add_info_to_file(input_file2, patient_clinics_dict)
    # add_index_to_file(input_file1)


if __name__ == '__main__':
    __main__()

