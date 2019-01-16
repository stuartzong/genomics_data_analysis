#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice

def make_patient_files_dict(bam_vcf_files):
    """ Dictionary holds all files: patient -> status -> file_identifier -> file_path  """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            patient = line['patient']
            status = line['status']
            lib = line['DNA_lib']
            DNA_tc = line['DNA_tumour_content']
            RNA_tc = line['RNA_tumour_content']
            DNA_single_vcf = line['DNA_single_vcf']
            paired_mpileup_vcf = line['paired_mpileup_vcf']
            mutseq_snv_vcf = line['mutseq_snv_vcf']
            DNA_bam = line['DNA_bam']
            strelka_snv_vcf = line['strelka_snv_vcf']
            strelka_indel_vcf = line['strelka_indel_vcf']
            other_vcf = line['other_vcf']
            cnv = line['cnv']
            RNA_bam = line['RNA_bam']
            RNA_single_vcf = line['RNA_single_vcf']
            tmp_dict = {'lib':lib,
                    'DNA_single_vcf':DNA_single_vcf,
                    'paired_mpileup_vcf': paired_mpileup_vcf,
                    'mutseq_snv_vcf': mutseq_snv_vcf,
                    'DNA_bam':DNA_bam,
                    'strelka_snv_vcf':strelka_snv_vcf,
                    'strelka_indel_vcf':strelka_indel_vcf,
                    'cnv':cnv,
                    'DNA_tc':DNA_tc,
                    'RNA_tc':RNA_tc,
                    'other_vcf':other_vcf,
                    'RNA_bam':RNA_bam,
                    'RNA_single_vcf':RNA_single_vcf}
            try:
                patient_files[patient][status].append(tmp_dict)
            except KeyError:
                if patient not in patient_files:
                    patient_files[patient] = {}
                if status not in patient_files[patient]:
                    patient_files[patient][status] = tmp_dict
        return patient_files




def make_bash_scripts(patient_files):
    cwdir = os.getcwd()
    somatic_dir = "/".join([cwdir, "somatic"])
    if not os.path.exists(somatic_dir):
        os.makedirs(somatic_dir)

    filter_script = "/home/szong/projects/development/variant/filterForSomatic.sh"
    for patient in patient_files:
        for status in patient_files[patient]:
            sub = '_'.join([patient, status])
            subdir = '/'.join([somatic_dir, sub])
            script = ".".join([subdir, 'sh'])
            #print dir
            if not os.path.exists(subdir):
                os.makedirs(subdir)
            vcfs = [patient_files[patient][status][i] for i in patient_files[patient][status] if ("vcf" in i and "other" not in i)]
            #print vcfs        
            with open(script, 'wb') as writer:
                writer.write("#! /bin/bash\n")
                #writer.write("cd %s/%s\n" % (somatic_dir, subdir))
                writer.write("cd %s\n" % subdir)
                for vcf in vcfs:
                    if (vcf !="NA"):
                        #print vcf.split('/')[-1]
                        writer.write("bash %s %s\n" % (filter_script, vcf))
                #writer.write("cd /projects/trans_scratch/validations/workspace/szong/David_Kaplan/variants/somatic/\n")

def __main__():
    print "scripts starts at: %s\n" % datetime.datetime.now()
    parser = argparse.ArgumentParser(description='Generate filter for somatic variant bash scripts')
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    input_file = args['input_file']


    print "Generating patient_files dictionary!\n"
    patient_files = make_patient_files_dict(input_file)

    make_bash_scripts(patient_files)

if __name__ == '__main__':
    __main__()

