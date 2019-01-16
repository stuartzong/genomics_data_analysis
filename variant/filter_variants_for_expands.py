#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice
import ConfigParser
import numpy as np

def __main__():
    print "Expand setup script starts at: %s\n" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='Filter variants based on qulaity and somatic filters')
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    parser.add_argument('-v','--variants_infile', help='specify quality filtered variant summary file', required=True)
    args = vars(parser.parse_args())

    # input_file = "bam_vcf_cnv_path.txt" 
    infile = args['input_file']
    filtered_summary = args['variants_infile']
    print "infile: %s" % infile
    print "filtered_summary: %s" % filtered_summary

    # make patient files dictionary
    patient_files = make_patient_files_dict(infile)
    out_dict = group_patients(patient_files)
    patient_files_wn = out_dict[0]
    patient_files_non = out_dict[1]
    


    print "Patients with matched normal are:"
    pprint(patient_files_wn)
    
    print "Patients without matched normal are:"
    pprint(patient_files_non)


    print "patient_biopsies dict is:"
    #variants_pat_dict = make_variant_dict(filtered_summary)
    patient_biopsies = make_biopsies_dict(patient_files)

    print "Dictionary indicates how many spatial/temporal samples have calls"
    variant_patients = make_variant_patients_dict(filtered_summary)
    pprint(variant_patients)

    print "Calculating alt_ploidy and put all afs for a variant in an array so that standard deviation can be computed"
    out = make_expand_dicts(filtered_summary, patient_files_wn, variant_patients, patient_biopsies)
    patient_variants = out[0] 
    variant_afs = out[1]
    pprint(patient_variants)
    pprint(variant_afs)

    print "Identify variants with differential afs among tumour spatial/temporal samples" 
    diff_variants = select_diff_variants(variant_afs)
    pprint(diff_variants)

    print "Writing expand summary file"
    make_expand_summary(filtered_summary, diff_variants)

    print "Making snv files"
    make_snv_files(patient_variants, diff_variants)
    print "Making cnv files"
    copy_cnv_files(patient_files_wn)
    print "Making R scripts"
    make_R_scripts(patient_files_wn)


def make_biopsies_dict(patient_files):
    biopsies = dict()
    for patient in patient_files:
        statuses = []
        for status in patient_files[patient]:
            if (status.lower() != 'normal'):
                statuses.append(status)
        num = len(statuses)
        biopsies[patient] = num
    print biopsies
    return biopsies



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



def group_patients(patient_files):
    patient_files_wn = dict()
    patient_files_non = dict()
    for patient in patient_files:
        tmp_dict = patient_files[patient]
        if ("normal" in patient_files[patient]):
            patient_files_wn[patient] = tmp_dict
        else:
            patient_files_non[patient] = tmp_dict
    return [patient_files_wn, patient_files_non] 

def make_variant_patients_dict(infile):
    #so that only variants with AFs available for all tumor tissues can be identified
    variant_patients = dict()
    with open (infile, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            gene = line['gene']
            chr = line["chromosome"]
            pos = line["position"]
            ref = line["ref_base"]
            alt = line["alt_base"]
            pat = line['patient_ID']

            try:
                DNA_t_af = float(line["t_DNA_AF"])
            except:
                 DNA_t_af = line["t_DNA_AF"]
            try:
                 DNA_n_af = float(line["n_DNA_AF"])
            except KeyError:
                DNA_n_af = "na"
            #remove variant on sex chromosomes or t_af < n_af
            if chr.isdigit() and (DNA_t_af > DNA_n_af):
               variant = '_'.join([gene, chr, pos, ref, alt])
               try:
                   variant_patients[variant].append(pat)
               except KeyError:
                   variant_patients[variant] = [pat]

    pprint(variant_patients)
    return variant_patients

def select_diff_variants(variant_afs):
    diff_variants = []
    # remove variants with indifferntial afs among biopsies 
    for variant in variant_afs:
        lst = [float(i) for i in variant_afs[variant]]
        std = np.std(lst)
        print lst, std
        if (std > 0.05):
            diff_variants.append(variant)
    print diff_variants
    return diff_variants         


def make_expand_summary(filtered_summary, diff_variants):
    expand_summary = ".".join([filtered_summary, "expand.txt"])
    with open (expand_summary,  'wb') as fh1:
        writer = csv.writer( fh1, delimiter='\t' )
        with open (filtered_summary, 'r')  as fh:
             records = csv.DictReader(fh,  delimiter='\t')
             headers = records.fieldnames
             writer.writerow(headers) 
             for line in records:
                 gene = line['gene']
                 chr = line["chromosome"]
                 pos = line["position"]
                 ref = line["ref_base"]
                 alt = line["alt_base"]
                 variant = '_'.join([gene, chr, pos, ref, alt])
                 if (variant in diff_variants): 
                     content = [line[i] for i in headers]
                     writer.writerow(content)




def make_expand_dicts(filtered_summary, patient_files_wn, variant_patients_dict, patient_biopsies):
    patient_variants = dict()
    variant_afs = {}
    # header = ["mutation_id", "ref_counts", "var_counts","normal_cn", "minor_cn", "major_cn", "variant_freq"]
    with open (filtered_summary, 'r')  as fh:
         records = csv.DictReader(fh,  delimiter='\t')
         for line in records:
             gene = line['gene']
             chr = line["chromosome"]
             pos = line["position"]
             ref = line["ref_base"]
             alt = line["alt_base"]
             DNA_n_altC = int(line["n_DNA_AltC"])
             DNA_n_af = float(line["n_DNA_AF"])
             DNA_t_refC = int(line["t_DNA_RefC"])
             DNA_t_altC = int(line["t_DNA_AltC"])
             DNA_t_af = float(line["t_DNA_AF"])
             adj_DNA_t_af = float(line["adj_t_DNA_AF"])
             pat = line["patient_ID"].split('_')
             patient = pat[0]
             status = '_'.join(pat[1:])
             variant = '_'.join([gene, chr, pos, ref, alt])
             if chr.isdigit() and (DNA_n_af < DNA_t_af): # tumour af has to be > n_af
                 if (len(variant_patients_dict[variant]) == patient_biopsies[patient]):
                     #make variant_afs so that standard deviation can be calculated for adj_DNA_t_afs
                     try:
                         variant_afs[variant].append(str(adj_DNA_t_af))
                     except KeyError:
                         variant_afs[variant] = [str(adj_DNA_t_af)]

                     #calculate alt_ploidy
                     if (DNA_n_af < 0.03 or DNA_n_altC < 2):
                         alt_ploidy = 0 # somatic 
                     else:
                         alt_ploidy = 1 # germline
                     variant = ":".join([chr, pos,  str(DNA_t_refC), str(DNA_t_altC), str(DNA_t_af), str(alt_ploidy), gene, ref, alt])
                     try:
                         patient_variants[patient][status].append( variant )
                     except KeyError:
                         if (patient not in patient_variants):
                             patient_variants[patient] = {}
                         if status not in patient_variants[patient]:
                             patient_variants[patient][status] = [variant]
    #print patient_variants, variant_afs
    return [patient_variants, variant_afs]




def make_snv_files(patient_variants, diff_variants): 
    snv_header = ["chr", "startpos", "T_refC", "T_altC","AF_Tumor", "PN_B", "gene"]
    for patient in patient_variants:
        afs = []
        for status in patient_variants[patient]:
            snv_file = "".join([patient, '_', status, ".expands.snv" ])
            with open (snv_file,  'wb') as fh:
                writer = csv.writer( fh, delimiter='\t' )
                writer.writerow( snv_header )
                for mutation in list(set(patient_variants[patient][status])):
                    sl =  mutation.split(":")
                    gene = sl[-3]
                    chr = sl[0]
                    pos = sl[1]
                    ref = sl[-2]
                    alt = sl[-1]
                    variant = '_'.join([gene, chr, pos, ref, alt])
                    #remove af no significant difference among various biopsies
                    if (variant in diff_variants):
                        writer.writerow(sl[:7])

def make_R_scripts(patient_files_wn):     
    #print "Generating expands R scripts!\n"
    wkdir = os.getcwd()
    wkdir = wkdir + "/"
    for patient in patient_files_wn:
        for status in patient_files_wn[patient]:
            if ("normal" not in status.lower()):
                blah = "_".join([patient, status])
                cn_file = "".join([wkdir, blah, ".expands.cn"])
                snv_file = "".join([wkdir, blah, ".expands.snv" ])
                pdf = "_".join([blah, 'expands'])
                dir ="".join([wkdir, "/", blah])
                r_script = ".".join([ blah, "r" ])
                with open ( r_script,  'wb' ) as info_fh:
                    info_fh.write( "".join(["setwd(\"",wkdir,"\")\n"]))
                    info_fh.write( "".join(["dir.create(\"", blah, "\", showWarnings = TRUE, recursive = FALSE, mode =""\"0777\")\n"]))
                    info_fh.write( "".join(["setwd(\"",dir,"\")\n"]))
                    info_fh.write( "library(expands)\n")
                    info_fh.write( "\n")
                    info_fh.write( "".join(["runExPANdS(\"", snv_file, "\",\"", cn_file, "\",maxScore=2.5, max_PM=6, min_CellFreq=0.1, precision=NA,plotF=2,snvF=\"", pdf, "\",maxN=8000,region=NA)\n"]) )
                    info_fh.write( "".join(["setwd(\"",wkdir,"\")\n"]))
 

def copy_cnv_files(patient_files_wn):
    #print "Copying cnv files!"
    cn_header = ["chr", "startpos", "endpos", "CN_Estimate"]
    wkdir = os.getcwd()
    wkdir = wkdir + "/"
    for patient in patient_files_wn:
        for status in patient_files_wn[patient]:
            if ("normal" not in status.lower()):
                blah = "_".join([patient, status])
                cn_file = "".join([wkdir, blah, ".expands.cn"])

                #copy cn file and add header, cn needs to be absolute copy number
                cnv = patient_files_wn[patient][status]["cnv"]
                shutil.copyfile(cnv, cn_file)
                for line in fileinput.input(cn_file, inplace=True):
                    if fileinput.isfirstline():
                        print '\t'.join(cn_header)
                    print line,



if __name__ == '__main__':
    __main__()

