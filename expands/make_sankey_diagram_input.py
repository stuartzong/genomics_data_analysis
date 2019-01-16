import os
import subprocess
import re
import sys
import glob
#import argparse
import csv
from collections import defaultdict
from pprint import pprint
import fileinput
import shutil

def make_patient_vcf_file_dict(bam_vcf_files):
    #patient>status>files
    cn_header = ["chr", "startpos", "endpos", "CN_Estimate"]

    patient_status_files = dict()
    with open (bam_vcf_files, 'r') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        tmp_dict=dict()
        for line in records:
            patient = "-".join(line['patient'].split("-")[0:3])
            status = line['status']
            lib = line['lib']
            tc = line['tumour_content']
            single_vcf = "NA"
            paired_vcf = line['paired_vcf']
            bam = line['bam_path']
            strelka_vcf = line['strelka_vcf']
            cnv = line['cnv']
            value=[lib, single_vcf, paired_vcf, strelka_vcf, bam, tc, cnv]
            #copy cn file and add header, cn needs to be absolute copy number
            if ("normal" not in status):
                blah = "_".join([patient, status])
                cn_file = ".".join([blah, "expands.cn"])
                shutil.copyfile(cnv, cn_file)
                for line in fileinput.input(cn_file, inplace=True):
                    if fileinput.isfirstline():
                        print '\t'.join(cn_header)
                    print line,
            
            try:
                #print patient, status
                patient_status_files[patient][status] = value
            except KeyError:
                #print "key error!"
                if patient not in patient_status_files:
                    patient_status_files[patient] = {}
                if status not in patient_status_files[patient]:
                    patient_status_files[patient][status] = value
        print patient_status_files
        return patient_status_files

def make_sankey_input_file( patient_files ):
    for patient in patient_files:
        #sps = sub populations
        primary_d = dict()
        refractory_d = dict()
        primary_snvd = dict()
        refractory_snvd = dict()
        wkdir = "/projects/trans_scratch/validations/workspace/szong/IF-AML/genome/expands/test/" 
        sankey_header = ["primary_sp", "mutation", "refractory_sp", "value" ]
        if (len(patient_files[patient]) == 3):
             all_mutations = []
             primary_sps = "".join([wkdir, patient, "_primary.expands.sps"])
             primary_snv= "".join([wkdir, patient, "_primary.expands.snv"])
             refractory_sps = "".join([wkdir, patient, "_malignant.expands.sps" ])
             refractory_snv = "".join([wkdir, patient, "_malignant.expands.snv" ])
             with open ( primary_sps, 'r' ) as fh1:
                 recd1 = csv.DictReader(fh1,  delimiter='\t')
                 primary_sps = []
                 for line in recd1:
                     chr = line["chr"]
                     startpos = line["startpos"]
                     sp = line["SP"]
                     mutation = "_".join([chr, startpos]) 
                     # NA means EXPANDS failed to compute cellular frequency for this SNV
                     if ( "NA" not in sp):
                         primary_sps.append(float(sp))
                         all_mutations.append(mutation)
                         try:
                             primary_d[mutation].append( sp )
                         except KeyError:
                             primary_d[mutation] = [sp]
             pprint(primary_d)
             max_primary_sp = max(primary_sps)
             print "max_primary_sp is:", max_primary_sp
            
             with open ( primary_snv, 'r' ) as fh3:
                 recd3 = csv.DictReader(fh3,  delimiter='\t')
                 for line in recd3:
                     chr = line["chr"]
                     startpos = line["startpos"]
                     gene = line["gene"]
                     mutation = "_".join([chr, startpos])
                     try:
                         primary_snvd[mutation].append( gene )
                     except KeyError:
                         primary_snvd[mutation] = [gene]
             pprint(primary_snvd)
 


             with open ( refractory_snv, 'r' ) as fh4:
                 recd4 = csv.DictReader(fh4,  delimiter='\t')
                 for line in recd4:
                     chr = line["chr"]
                     startpos = line["startpos"]
                     gene = line["gene"]
                     mutation = "_".join([chr, startpos])
                     try:
                         refractory_snvd[mutation].append( gene )
                     except KeyError:
                         refractory_snvd[mutation] = [gene]
             print "xxxxxxxxxxxxxxxx", patient
             pprint(refractory_snvd)



             with open ( refractory_sps, 'r' ) as fh2:
                 recd2 = csv.DictReader(fh2,  delimiter='\t')
                 refractory_sps = []
                 for line in recd2:
                     chr = line["chr"]
                     startpos = line["startpos"]
                     sp = line["SP"]
                     mutation = "_".join([chr, startpos])
                     # NA means EXPANDS failed to compute cellular frequency for this SNV
                     if ( "NA" not in sp):
                         refractory_sps.append(float(sp))
                         all_mutations.append(mutation)
                         try:
                             refractory_d[mutation].append( sp )
                         except KeyError:
                             refractory_d[mutation] = [sp]
             pprint(refractory_d)
             all_mutations = list(set(all_mutations))
             print "all_mutations are:", all_mutations 
             max_refractory_sp = max(refractory_sps)
             print "max_refractory_sp is:", max_refractory_sp

             #adjust cellular frequency to 100% tumour content, write to sankey input file
             sankey_file = "".join([wkdir, patient, ".sankey" ])
             with open ( sankey_file,  'wb' ) as fh3:
                 writer = csv.writer( fh3, delimiter='\t' )
                 writer.writerow( sankey_header )
                 for mutation in all_mutations:
                     try:
                         primary_sp = primary_d[mutation][0]
                         if ("NA" not in primary_sp):
                             primary_sp = "{:.3f}".format(float(primary_sp)/max_primary_sp)
                             primary_sp = "_".join(["primary",str(primary_sp)])
                     except KeyError:
                         primary_sp = "primary_NA"
                     try:
                         refractory_sp = refractory_d[mutation][0]
                         if ("NA" not in refractory_sp):
                             refractory_sp = "{:.3f}".format(float(refractory_sp)/max_refractory_sp)
                             refractory_sp = "_".join(["refractory",str(refractory_sp)] )
                     except KeyError:
                         refractory_sp = "refractory_NA"
                     try:
                         geneX = primary_snvd[mutation][0]
                     except KeyError:
                         geneX = refractory_snvd[mutation][0]
                     mutation = "_".join([geneX, mutation])
                     sankey_width = 0.075
                     if ("MODIFIER" not in mutation):
                         sankey_width = 0.2 
                     writer.writerow([primary_sp, mutation, refractory_sp, sankey_width] )
     

def __main__():

    print "\nGenerating patient_files dictionary!"
    input_files = "bam_vcf_cnv_path.txt"
    patient_files = make_patient_vcf_file_dict(input_files)
    pprint(patient_files)

    #variant_summary = "SNV_summary_with_normal.txt.filtered"
    #variant_summary = "SNV_summary_with_normal.txt.all"
    make_sankey_input_file(patient_files)

if __name__ == '__main__':
    __main__()


