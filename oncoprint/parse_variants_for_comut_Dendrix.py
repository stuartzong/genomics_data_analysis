#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice
import ConfigParser

def __main__():
    print "script starts at: %s\n" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='Reformat filtered summary for comut plot!')
    parser.add_argument('-i1','--input_file1', help='specify SNV summary file', required=True)
    parser.add_argument('-i2','--input_file2', help='specify INDEL summary file', required=True)
    args = vars(parser.parse_args())
    variant_summary = args['input_file1']
    indel_summary = args['input_file2']

    gene_list = "ordered_genes.txt"
    #patient_list = "HIV_Cervical_patients.txt"
    patient_list = "HIV_Cervical_101_patients.txt.hiv.sorted"
    #patient_list = "HIV_Cervical_101_patients.txt.histology.sorted"
    # variant_input_file = "bam_vcf_cnv_path.txt" 
    comut_summary = ".".join([variant_summary, "comut", "hiv", "sorted" ])

    genes = make_list(gene_list)
    patients = make_list(patient_list)

    format_variants(variant_summary, indel_summary, comut_summary, genes, patients)
    # pat_vars = variants_4_Dendrix(variant_summary, genes)
    # Dendrix_file = 'Dendrix_input_file.txt'
    # write_Dendrix_file(pat_vars, Dendrix_file)

def make_list(infile):
    print "infile is: ", infile
    items = []
    with open(infile, 'r') as fh:
        for line in fh:
            item = line.split()[0]
            items.append(item)
        #items = list(set(items))
        print items
    return items


def variants_4_Dendrix(variant_summary, genes):
    pat_vars = dict()
    with open(variant_summary, 'r') as fh:
        records = csv.DictReader(fh, delimiter='\t')
        headers = records.fieldnames
        for line in records:
            patient = line["patient_ID"].split('_')[0]
            gene = line['gene']
            if gene in genes:
                try:
                    pat_vars[patient].append(gene)
                except KeyError:
                    pat_vars[patient] = [gene]

    for pat in pat_vars:
        genes = pat_vars[pat]
        pat_vars[pat] = list(set(genes))
    pprint(pat_vars)
    return pat_vars

def write_Dendrix_file(pat_vars, outfile):
    with open(outfile, 'wb') as opf:
        writer = csv.writer(opf, delimiter='\t')
        for pat in pat_vars:
            genes = pat_vars[pat]
            content = [pat] + genes
            writer.writerow(content)


def format_variants(variant_summary, indel_summary, comut_summary, genes, patients):

    """ d = gene_variant_patients dictionary """
    d = dict()
    with open(comut_summary,  'wb') as fh:
            headers = ['gene', 'sample', 'mutations']
            writer = csv.writer( fh, delimiter='\t' )
            writer.writerow(headers)
            print "The variant summary file is: %s.\n" % variant_summary
            with open(variant_summary, 'r') as handle:
                records = csv.DictReader(handle,  delimiter='\t')
                headers = records.fieldnames
                for line in records:
                    gene = line['gene']
                    chr = line["chromosome"]
                    pos = line["position"]
                    ref = line["ref_base"]
                    alt = line["alt_base"]
                    patient = line["patient_ID"].split('_')[0]
                    # fix:some snpeff noted as lowercase, which gives trouble for comut plot
                    mutation= re.split('\(|', line["snpeff_details"])[0].upper()
                    # be aware that there are annotations as 
                    # MISSENSE_VARIANT+SPLICE_REGION_VARIANT
                    # SPLICE_ACCEPTOR_VARIANT+INTRON_VARIANT
                    if "+" in mutation:
                        print "xxxxxxxx", mutation


                    try:
                        d[gene][patient].append(mutation)
                    except KeyError:
                        if gene not in d:
                            d[gene] = {}
                        if patient not in d[gene]:
                            d[gene][patient] = [mutation]



            # add indels into the comut summary 
            with open(indel_summary, 'r') as indel_handle:
                records = csv.DictReader(indel_handle, delimiter='\t')
                headers = records.fieldnames
                for line in records:
                    gene = line['gene']
                    ref = line['ref_base']
                    alt = line['alt_base']
                    patient = line['patient_ID'].split('_')[0]
                    mutation= re.split('\(|', line["snpeff_details"])[0].upper()
                    try:
                        d[gene][patient].append(mutation)
                    except KeyError:
                        if gene not in d:
                            d[gene] = {}
                        if patient not in d[gene]:
                            d[gene][patient] = [mutation]





            #for gene in d:
            #    if (gene in genes):
            for gene in genes:
                if gene in d: 
                    if (len(d[gene]) > 0):
                        for patient in patients:
                            # for  patient in d[gene]:
                            if (patient in d[gene]):
                                mut = list(set(d[gene][patient]))
                                if (int(len(mut)) > 1):
                                    mutation = "Multiple"
                                else:
                                    mutation = mut[0]
                            else:
                                mutation = "NA"
                                mut = ["NA"]
                            print gene, "\t",patient, "\t", mutation,"\t", len(mut)
                            writer.writerow([gene, patient, mutation])

                            
if __name__ == '__main__':
    __main__()

