#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice

# test comments
def make_patient_files_dict(bam_vcf_files):
    """ Dictionary holds all files: patient -> status -> file_identifier -> file_path  """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            patient = line['patient']
            status = line['status']
            DNA_bbt = line['bbt_genome']
            RNA_bbt = line['bbt_transcriptome']
            HIV_status = line['HIV_status']
            tmp_dict = {'HIV_status':HIV_status,
                    'DNA_bbt': DNA_bbt,
                    'RNA_bbt': RNA_bbt}
            try:
                patient_files[patient][status].append(tmp_dict)
            except KeyError:
                if patient not in patient_files:
                    patient_files[patient] = {}
                if status not in patient_files[patient]:
                    patient_files[patient][status] = tmp_dict
        return patient_files



def is_other_readable(file_path):
    m = os.stat(file_path).st_mode
    #otherExec  = bool(m & 0001)
    #otherWrite = bool(m & 0002)
    other_read_status  = bool(m & 0004)
    #group_read_status = bool(m & stat.S_IRGRP)
    #return group_read_status
    return other_read_status


def exist(file_path):
    exist_status = os.path.exists(file_path)
    return exist_status
 
def check_files(files):
    missing_files = []
    if (len(files) == 0):
       print "No files to check!"
    for file in files:
        short_name = file.split("/")[-1]
        exist_status = False
        read_status = False
        exist_status = exist(file)
        if (exist_status):
            read_status = is_other_readable(file)
        if (not(exist_status and read_status)):
            missing_files.append(file)
        else:
            print "%s: OK" % short_name 
    if missing_files:
        print "ERROR: The following files are either missing or no read permission! Check if pipeline anlyais has completed!"
        for file in missing_files:
            print "Missing %s\n" % file
            sys.exit()



def parse_bbt_file(infile):
    microbe_counts = dict()
    with open (infile, 'r') as fh:
        for line in fh:
            sl = line.strip().split()
            reads = sl[0]
            microbe = sl[1]
            microbe_counts[microbe] = reads
    #pprint(microbe_counts)
    return microbe_counts





def check_file_permission(patient_files):
    DNA_bbts = []
    RNA_bbts = []
    for patient in patient_files:
        for status in patient_files[patient]:
            DNA_bbt = patient_files[patient][status]['DNA_bbt']
            RNA_bbt = patient_files[patient][status]['RNA_bbt']

            if (DNA_bbt != "NA"):
                DNA_bbts.append(DNA_bbt)
            if (RNA_bbt != "NA"):
                RNA_bbts.append(RNA_bbt)

    print "Checking DNA bbt file permissions!"
    check_files(DNA_bbts)

    print "Checking RNA bbt file permissions!"
    check_files(RNA_bbts)



def make_microbe_list(infile):
    microbes = []
    with open (infile, 'r')  as fh:
         records = csv.DictReader(fh,  delimiter='\t')
         for line in records:
             microbe = line['filter_id'].split('.')[0]
             microbes.append(microbe)
    print microbes
    return microbes 

def __main__():
    print "BBT summary scripts starts at: %s\n" % datetime.datetime.now()     

    parser = argparse.ArgumentParser(description='generate BBT result matrix')
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    args = vars(parser.parse_args())

    # input file contains the full path to *.adjusted.normalized.txt files
    input_file = args['input_file']

    print "BBT input file is:\n%s\n" % (input_file)

    print "Generating patient_files dictionary!\n" 
    patient_files = make_patient_files_dict(input_file)

    print "Checking bbt file permissions!"
    check_file_permission(patient_files)

    print "Making bbt filters list!"
    filters_file = "filters.txt"
    microbes = make_microbe_list(filters_file)


    print "Parsing bbt files"
    outfile = 'BBT_matrix.txt'
    with open (outfile, 'wb') as fh:
        writer = csv.writer(fh, delimiter='\t')
        headers = ['patient', 'tissue_status', 'HIV_status', 'data_type'] + microbes 
        writer.writerow(headers)
        for patient in patient_files:
            DNA_microbes = {}
            RNA_microbes = {}
            for status in patient_files[patient]:
                DNA_bbt = patient_files[patient][status]['DNA_bbt']
                RNA_bbt = patient_files[patient][status]['RNA_bbt']
                HIV_status = patient_files[patient][status]['HIV_status']
                if (DNA_bbt != 'NA'):
                    DNA_microbes = parse_bbt_file(DNA_bbt)
                if (RNA_bbt != 'NA'):
                    RNA_microbes = parse_bbt_file(RNA_bbt)
                # print out microbes results
                DNA_count = [ format(float(DNA_microbes[microbe].strip('%'))/100, '.10f') if microbe in DNA_microbes else 0 for microbe in microbes ]
                RNA_count = [ format(float(RNA_microbes[microbe].strip('%'))/100, '.10f') if microbe in RNA_microbes else 0 for microbe in microbes ]
                print patient, status, HIV_status, DNA_count, RNA_count
                DNA_content = [patient, status, HIV_status, 'genome'] + DNA_count
                RNA_content = [patient, status, HIV_status, 'transcriptome'] + RNA_count
                writer.writerow(DNA_content)
                writer.writerow(RNA_content)

    print "Summarization scripts finished on: %s\n" % datetime.datetime.now()    

if __name__ == '__main__':
    __main__()

