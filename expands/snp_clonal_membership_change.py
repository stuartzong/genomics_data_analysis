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

def make_sps_dict( patient_files ):
    #sps_d = cancer cell subpopulations dictionary
    #cf = cellular frequency
    sps_d = dict()
    snp_d = dict()
    snp_cfchange_d = dict()
    for patient in patient_files:
        wkdir = os.getcwd()
        wkdir = wkdir + "/"
        #wkdir = "/projects/trans_scratch/validations/workspace/szong/IF-AML/genome/expands/test/" 
        if (len(patient_files[patient]) == 3):
             all_mutations = []
             sankey = "".join([wkdir, patient, ".sankey"])
             with open ( sankey, 'r' ) as fh1:
                 recd1 = csv.DictReader(fh1,  delimiter='\t')
                 primary_sps = []
                 #l2h = "l2h"
                 #h2l = "h2l"
                 #p2a = "p2a"
                 #a2p = "a2p"
                 #d2d = "d2d"
                 for line in recd1:
                     primary_sp = line["primary_sp"]
                     refractory_sp = line["refractory_sp"]
                     mutation = line["mutation"]
                     if ( "NA" not in primary_sp) and ( "NA" not in refractory_sp):
                         primary_cf = float(primary_sp.split("_")[1])
                         refractory_cf = float(refractory_sp.split("_")[1])
                         if (primary_cf >= 0.99 and refractory_cf >= 0.99):
                             key = "d2d"
                         elif (primary_cf > refractory_cf):
                             key = "h2l"
                         elif (primary_cf < refractory_cf):
                             key = "l2h"
                         # check if cf increase or decrease
                         if (refractory_cf >= primary_cf):
                             keyb = "resistant"
                         else:
                             keyb = "sensitive"
                     elif (  "NA" in primary_sp ):
                         key = "a2p"
                         keyb= "resistant"
                         primary_cf = 0
                         refractory_cf = float(refractory_sp.split("_")[1])
                     elif (  "NA" in refractory_sp ):
                         key = "p2a"
                         keyb = "sensitive"
                         refractory_cf = 0
                         primary_cf = float(primary_sp.split("_")[1])
                     else:
                         print "ERROR! invalid sp entry!"
                         primary_cf = "weird"
                         refractory_cf = "weird_too"
                     # dictionary: patient >> cf change type >> snps
                     try:
                         sps_d[patient][key].append( mutation )
                     except KeyError:
                         if patient not in sps_d:
                             sps_d[patient] = {}
                         if key not in sps_d[patient]:
                             sps_d[patient][key] = [mutation]
                     # dictionary: snp >> cf change type >> patient
                     try:
                         snp_d[mutation][key].append( patient )
                     except KeyError:
                         if mutation not in snp_d:
                             snp_d[mutation] = {}
                         if key not in snp_d[mutation]:
                             snp_d[mutation][key] = [patient]
 
                     # dictionary for resistant and sensitive mutations
                     details = [patient, primary_cf, refractory_cf]
                     try:
                         snp_cfchange_d[mutation][keyb].append( details )
                     except KeyError:
                         if mutation not in snp_cfchange_d:
                             snp_cfchange_d[mutation] = {}
                         if key not in snp_cfchange_d[mutation]:
                             snp_cfchange_d[mutation][keyb] = [ details ]
                           

    #pprint(sps_d)
    pprint(snp_d)
 
    #print out snp_d
    out_file = "snp_membership_change.txt"
    header = ["mutation", "numberTypesOfCellularFrequencyChange", "TypeOfCellularFrequencyChange", "numberOfRefractoryPatients", "RefractoryPatientID", "primaryCellularFrequency", "refractoryCellularFrequency"]
    with open ( out_file,  'wb' ) as fh3:
        writer = csv.writer( fh3, delimiter='\t' )
        writer.writerow( header )

        for mutation in snp_d:
            num_change = str(len(snp_d[mutation].keys()))
            for key in snp_d[mutation]:
                patients = snp_d[mutation][key]
                num_pat = str(len(patients))
                for patient in patients:
                    print "\t".join([mutation, num_change, key, num_pat, patient ]) 
                    writer.writerow([mutation, num_change, key, num_pat, patient] )
 

   #print out snp_cfchange_d
    out_file = "snp_cf_change.txt"
    #header = ["mutation", "numberTypesOfCellularFrequencyChange", "TypeOfMutation", "numberOfPatients", "patient"]
    with open ( out_file,  'wb' ) as fh4:
        writer = csv.writer( fh4, delimiter='\t' )
        writer.writerow( header )

        for mutation in snp_cfchange_d:
            num_change = str(len(snp_cfchange_d[mutation].keys()))
            for keyb in snp_cfchange_d[mutation]:
                patients = snp_cfchange_d[mutation][keyb]
                num_pat = str(len(patients))
                for patient in patients:
                    #print "\t".join([mutation, num_change, keyb, num_pat, patient[0], patient[1], patient[2] ])
                    writer.writerow([mutation, num_change, keyb, num_pat, patient[0], str(patient[1]), str(patient[2])] )



def __main__():

    print "\nGenerating patient_files dictionary!"
    input_files = "bam_vcf_cnv_path.txt"
    patient_files = make_patient_vcf_file_dict(input_files)
    pprint(patient_files)

    make_sps_dict(patient_files)

if __name__ == '__main__':
    __main__()


