#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice

# test comments
"""
testing script: expands pipeline (v1.0) for production
To run this script, type the following command:
python  variant_summarization_clonality_analysis.py -dt wgs or spc > summary_log_file.txt
normal: matched normal samples must be specified as "normal", case sensitive
"""
def make_patient_files_dict(patient_bam_file):
    patient_bam = dict()
    with open(patient_bam_file, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            patient = line['patient']
            bam = line['bam']
            patient_bam[patient] = bam
        return patient_bam

def is_group_readable(file_path):
    m = os.stat(file_path).st_mode
    group_read_status = bool(m & stat.S_IRGRP)
    return group_read_status


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
            read_status = is_group_readable(file)
        if (not(exist_status and read_status)):
            missing_files.append(file)
        else:
            print "%s: OK" % short_name 
    if missing_files:
        print "ERROR: The following files are either missing or no read permission! Check if STRELKA has completed!"
        for file in missing_files:
            print "Missing %s\n" % file
            sys.exit()

def make_cp_scripts(patient_bam):
    reference_genome= "/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa"
    samtools_path = "/home/rcorbett/aligners/samtools/samtools-0.1.17/samtools"
    bash_files = []
    stamps = []
    for patient in patient_bam:
        bam = patient_bam[patient]
        cp_bam = ".".join([patient, 'bam'])
        bash_file = ".".join([patient, "cpbam.sh"])
        bash_files.append(bash_file)
        stamp = "_".join([patient, 'success'])
        stamps.append(stamp)
        with open(bash_file,  'wb') as writer:
             writer.write("#! /bin/bash\n")
             writer.write("#$ -S /bin/bash\n")
             writer.write("#$ -N copybam\n")
             writer.write("#$ -q transabyss.q\n")
             writer.write("#$ -l mem_token=2G,mem_free=2G,h_vmem=2G\n")
             writer.write("#$ -V\n")
             writer.write("cp %s %s\n" % (bam, cp_bam))
             writer.write(" ".join([samtools_path, "index", cp_bam,"\n"]))
             writer.write("if [ $? -eq 0 ]\n")
             writer.write("    then\n")
             writer.write("    touch %s\n" % (stamp)) 
             writer.write("else\n")
             writer.write("    echo \"ERROR: cp and/or samtools index did not finish correctly!\" >>\"summary_log_file.txt\"\n")
             writer.write("    echo \"Failed patient is: %s \" >> \"summary_log_file.txt\"\n" % (patient)) 
             writer.write("    exit\n")
             writer.write("fi\n")
    return [bash_files, stamps]



def qsub_scripts(scripts):
    """ qsub scripts """
    wkdir = os.getcwd()
    for script in scripts:
        #p = subprocess.Popen('qsub %s' % script,  shell=True, stdout=subprocess.PIPE)
        p = subprocess.Popen('ssh tachpc \"cd %s;  qsub %s\"' % (wkdir, script),  shell=True, stdout=subprocess.PIPE)
        output,  err = p.communicate()



def detect_cluster_job_status(completeion_file_list):
    completed = False
    for file in completeion_file_list:
        if (os.path.exists(file)):
            completed = True
        else:
            completed = False
            break
    return completed

def detect_cluster_jobs(complete_stamps):
    """ detect if job on cluster finised """
    job_status = False
    print "Waiting for cluster jobs to finish!\n"
    while (not job_status):
        time.sleep(10)
        job_status = detect_cluster_job_status(complete_stamps)
    print "All cluster jobs finished? %s\n" % job_status


def check_file_permission(patient_files):
    bams = []
    for patient in patient_files:
            bam = patient_files[patient]
            if (bam != "NA"):
                bams.append(bam)
    print "Checking bam files permissions!"
    check_files(bams)



     
def __main__():
    print "Copying and indexing scripts starts at: %s\n" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='copy and index bam files')
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    patient_bam = args['input_file']


    patient_files = make_patient_files_dict(patient_bam)
    check_file_permission(patient_files)
    
    
    print "Generating cp scripts!"
    out = make_cp_scripts(patient_files) 
    cp_scripts = out[0]
    cp_complete_stamps = out[1]
    print cp_scripts

    qsub_scripts(cp_scripts)

    
    print "Detecting if cluster jobs finised!\n"
    detect_cluster_jobs(cp_complete_stamps)
            

if __name__ == '__main__':
    __main__()

