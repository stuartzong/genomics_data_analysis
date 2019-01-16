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
    parser.add_argument('-f','--quality_somatic_filters', help='specify quality and somatic filters', required=True)
    parser.add_argument('-p','--pairing', help='specify if sample paired with matched normal: paired or unpaired', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    variant_input_file = args['input_file']
    filters = args['quality_somatic_filters']
    pairing = args['pairing']
    filter_variants(variant_input_file, filters, pairing)

def filter_variants(variant_summary, filters, pairing):
    filtered_summary = ".".join([variant_summary, "filtered" ])
    somatic_summary = ".".join([filtered_summary, "somatic", "txt" ])

    # get filter values
    config = ConfigParser.SafeConfigParser()
    config.read(filters)


    # default quality filtering: 
    DNA_t_cov_thres = float(config.get('quality_filters', 'DNA_t_cov'))
    DNA_t_altC_thres = float(config.get('quality_filters', 'DNA_t_altC'))
    DNA_t_af_thres = float(config.get('quality_filters', 'DNA_t_af'))

    # OR 
    RNA_t_cov_thres = float(config.get('quality_filters', 'RNA_t_cov'))
    RNA_t_altC_thres = float(config.get('quality_filters', 'RNA_t_altC'))
    RNA_t_af_thres = float(config.get('quality_filters', 'RNA_t_af'))
    #misalignment at exon junctions
    RNA_t_percent_thres = float(config.get('quality_filters', 'RNA_t_altref_total_percent'))
    # third allele
    DNA_t_percent_thres = float(config.get('quality_filters', 'DNA_t_altref_total_percent'))


    # default somatic filters
    DNA_n_af_thres = float(config.get('somatic_filters','DNA_n_af'))
    DNA_n_altC_thres = float(config.get('somatic_filters', 'DNA_n_altC'))
    RNA_n_af_thres = float(config.get('somatic_filters', 'RNA_n_af'))
    RNA_n_altC_thres = float(config.get('somatic_filters', 'RNA_n_altC'))

    """ d = gene_variant_patients dictionary """
    d = dict()
    with open (filtered_summary,  'wb') as fh:
        with open (somatic_summary,  'wb') as fh2:
            writer = csv.writer( fh, delimiter='\t' )
            writer2 = csv.writer( fh2, delimiter='\t' )
            print "The variant summary file is: %s.\n" % variant_summary
            with open (variant_summary, 'r') as handle:
                 records = csv.DictReader(handle,  delimiter='\t')
                 headers = records.fieldnames
                 writer.writerow(headers)
                 writer2.writerow(headers)
                 for line in records:
                     gene = line['gene']
                     chr = line["chromosome"]
                     pos = line["position"]
                     ref = line["ref_base"]
                     alt = line["alt_base"]
                     try:
                         DNA_t_cov = int(line["t_DNA_cov"])
                     except:
                         DNA_t_cov = line["t_DNA_cov"]
                     try:
                         DNA_t_refC = int(line["t_DNA_RefC"])
                     except:
                         DNA_t_refC = line["t_DNA_RefC"]
                     try:
                         DNA_t_altC = int(line["t_DNA_AltC"])
                     except:
                         DNA_t_altC = line["t_DNA_AltC"]
                     try:
                         DNA_t_af = float(line["t_DNA_AF"])
                     except:
                         DNA_t_af = line["t_DNA_AF"]
                     try:
                         RNA_t_cov = int(line["t_RNA_cov"])
                     except:
                         RNA_t_cov = line["t_RNA_cov"]
                     try:
                         RNA_t_refC = int(line["t_RNA_RefC"])
                     except:
                         RNA_t_refC = line["t_RNA_RefC"]
                     try:
                         RNA_t_altC = int(line["t_RNA_AltC"])
                     except:
                         RNA_t_altC = line["t_RNA_AltC"]
                     try:
                         RNA_t_af = float(line["t_RNA_AF"])
                     except:
                         RNA_t_af = line["t_RNA_AF"]
                     RNA_t_altref_total = RNA_t_refC + RNA_t_altC
                     DNA_t_altref_total = DNA_t_refC + DNA_t_altC

                     try:  
                         DNA_n_cov = int(line["n_DNA_cov"])
                         DNA_n_refC = int(line["n_DNA_RefC"])
                         DNA_n_altC = int(line["n_DNA_AltC"])
                         DNA_n_af = float(line["n_DNA_AF"])
                     except KeyError:
                         DNA_n_cov = "na"
                         DNA_n_refC = "na"
                         DNA_n_altC = "na"
                         DNA_n_af = "na"
                     try:
                         RNA_n_cov = int(line["n_RNA_cov"])
                         RNA_n_refC = int(line["n_RNA_RefC"])
                         RNA_n_altC = int(line["n_RNA_AltC"])
                         RNA_n_af = float(line["n_RNA_AF"])
                     except:
                         RNA_n_cov = line["n_RNA_cov"]
                         RNA_n_refC = line["n_RNA_RefC"]
                         RNA_n_altC = line["n_RNA_AltC"]
                         RNA_n_af = line["n_RNA_AF"]
                         #RNA_n_cov = "na"
                         #RNA_n_refC = "na"
                         #RNA_n_altC = "na" 
                         #RNA_n_af = "na"
                           
                     content = [line[i] for i in headers]
                     #print type(DNA_t_cov_thres), DNA_t_cov_thres
                     #quality filtering
                     if (DNA_t_cov == 'na'):
                        print "Only transcriptome is sequenced for this tumor! "
                        if (RNA_t_cov >= RNA_t_cov_thres and RNA_t_altC >= RNA_t_altC_thres and RNA_t_af >= RNA_t_af_thres and RNA_t_altref_total >= RNA_t_percent_thres*RNA_t_cov):
                             writer.writerow(content)
                             print "aaaaaaaaaaa"
                             if ((RNA_n_af <= RNA_n_af_thres) or (RNA_n_altC <= RNA_n_altC_thres)):
                                 writer2.writerow(content)
                     elif (RNA_t_cov == 'na'):
                         print "Only genome is sequenced for this tumor! "
                         if (DNA_t_cov >= DNA_t_cov_thres and DNA_t_altC >= DNA_t_altC_thres and DNA_t_af >= DNA_t_af_thres and DNA_t_altref_total >= DNA_t_percent_thres*DNA_t_cov):
                             writer.writerow(content)
                             print "bbbbbbbb"
                             if ((DNA_n_af <= DNA_n_af_thres) or (DNA_n_altC <= DNA_n_altC_thres)):
                                 writer2.writerow(content)
 
                     else:
                         print "Both genome and transcriptome are sequenced for this tumor! "
                         if ((DNA_t_cov >= DNA_t_cov_thres and DNA_t_altC >= DNA_t_altC_thres and DNA_t_af >= DNA_t_af_thres and DNA_t_altref_total >= DNA_t_percent_thres*DNA_t_cov) or
                             (RNA_t_cov >= RNA_t_cov_thres and RNA_t_altC >= RNA_t_altC_thres and RNA_t_af >= RNA_t_af_thres and RNA_t_altref_total >= RNA_t_percent_thres*RNA_t_cov)):
                             writer.writerow(content)
                             
                             #if (pairing == "paired"):
                             # somatic filters
                             if (DNA_n_cov == 'na'):
                                 if ((RNA_n_af <= RNA_n_af_thres) or (RNA_n_altC <= RNA_n_altC_thres)):
                                     writer2.writerow(content)
                             elif (RNA_n_cov == 'na'):
                                 if ((DNA_n_af <= DNA_n_af_thres) or (DNA_n_altC <= DNA_n_altC_thres)):
                                     writer2.writerow(content)
                             else:
                                 if ((DNA_n_af <= DNA_n_af_thres) or (DNA_n_altC <= DNA_n_altC_thres)) and ((RNA_n_af <= RNA_n_af_thres) or (RNA_n_altC <= RNA_n_altC_thres)):
                                     writer2.writerow(content)
    return [filtered_summary, somatic_summary]


if __name__ == '__main__':
    __main__()

