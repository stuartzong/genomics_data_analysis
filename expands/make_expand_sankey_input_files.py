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

def __main__():
    print "Generating patient_files dictionary!\n"
    input_files = "NB_TIC_128153295_lines_for_BC_genome_normal_cn_unknown.csv"
    patient_files = make_patient_files_dict(input_files)
    #variant_summary = "high_moderate_SNV_summary_no_normal.txt.128153295.somatic.filtered"
    #make_expands_input_file(variant_summary, patient_files)
    make_sankey_input_file(patient_files)

def copy_cn_file(patient, status, cn_file):
    """ Copy cn file and add header, cn has to be absolute copy number! """
    cn_header = ["chr", "startpos", "endpos", "CN_Estimate"]
    pat_status = "_".join([patient, status])
    local_cn_file = ".".join([pat_status, "expands.cn"])
    shutil.copyfile(cn_file, local_cn_file)
    for line in fileinput.input(local_cn_file, inplace=True):
        if fileinput.isfirstline():
            print '\t'.join(cn_header)
        print line,


def make_patient_files_dict(bam_vcf_files):
    """
    Generate a dictionary: patient -> status -> 
    [lib, single_vcf, paired_vcf, strelka_snv_vcf, bam, tc, cnv, strelka_indel_vcf]
    """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            patient = line['patient']
            status = line['status']
            lib = line['DNA_lib']
            tc = line['tumour_content']
            single_vcf = line['single_vcf']
            paired_vcf = line['paired_mpileup_vcf']
            bam = line['bam_path']
            strelka_snv_vcf = line['strelka_snv_vcf']
            strelka_indel_vcf = line['strelka_indel_vcf']
            cn_file = line['cnv']
            value = [lib, single_vcf, paired_vcf, strelka_snv_vcf, bam, tc, cn_file, strelka_indel_vcf]
            if ("normal" not in status):
                copy_cn_file(patient, status, cn_file)

            try:
                patient_files[patient][status].append(value)
            except KeyError:
                if patient not in patient_files:
                    patient_files[patient] = {}
                if status not in patient_files[patient]:
                    patient_files[patient][status] = value
        pprint(patient_files)
        return patient_files

def make_expands_input_file( variant_summary, patient_files ):
    """ 
    d = gene_variant_patients dictionary 
    tc_d = tumor content dictionary
    """
    d = dict()
    tc_d = dict()
    snv_header = ["chr", "startpos", "T_refC", "T_altC","AF_Tumor", "PN_B", "gene"]
    # header = ["mutation_id", "ref_counts", "var_counts","normal_cn", "minor_cn", "major_cn", "variant_freq"]
    with open ( variant_summary, 'r' ) as handle:
         records = csv.DictReader(handle,  delimiter='\t')
         for line in records:
             gene = line["gene"]
             chr = line["chromosome"]
             pos = line["position"]
             ref = line["referenceBase"]
             alt = line["alternativeBase"]
             pat = line["patientID"]
             t_refC = line["tumourDNARefBaseCount"]
             t_altC = line["tumourDNAAltBaseCount"]
             t_tot =line["tumourDNASequencingCoverage"]
             t_af = line["tumourDNAAlleleFrequency"]
             if chr.isdigit() and (float(t_af) > 0.1): # tumour af has to be > n_af
                 ''' assume all variants are somatic
                 if ("not_in_cosmic64" not in cosm_id):
                     alt_ploidy = 0
                 elif ("novel_snp" in SNPID): 
                     alt_ploidy = 0 # somatic 
                 else:
                     alt_ploidy = 1 # germline
                 '''
                 alt_ploidy = 0 # somatic
                 variant = ":".join([chr, pos,  t_refC, t_altC, t_af, str(alt_ploidy), gene])
                 try:
                     d[pat].append( variant )
                 except KeyError:
                     #print "key error!"
                     d[pat] = [variant]
    print d

    print "Making expands tsv input files!\n"
    for pat in d:
        snv_file = "".join([pat, ".expands.snv" ])
        with open ( snv_file,  'wb' ) as fh:
            writer = csv.writer( fh, delimiter='\t' )
            writer.writerow( snv_header )
            for mutation in list(set(d[pat])):
                sl =  mutation.split(":")
                writer.writerow(sl)

    print "Generating expands R scripts!\n" 
    wkdir = os.getcwd() 
    wkdir = wkdir + "/"
    for patient in patient_files:
        for status in patient_files[patient]:
            if ("normal" not in status):
                blah = "_".join([patient, status])
                cn_file = "".join([wkdir, blah, ".expands.cn"])
                snv_file = "".join([wkdir, blah, ".expands.snv" ])
                dir ="".join([wkdir, "/", blah])
                r_script = ".".join([ blah, "r" ]) 
                with open ( r_script,  'wb' ) as info_fh:
                    info_fh.write( "".join(["setwd(\"",wkdir,"\")\n"]))
                    info_fh.write( "".join(["dir.create(\"", blah, "\", showWarnings = TRUE, recursive = FALSE, mode =""\"0777\")\n"]))
                    info_fh.write( "".join(["setwd(\"",dir,"\")\n"]))
                    info_fh.write( "library(expands)\n")
                    info_fh.write( "\n")
                    info_fh.write( "".join(["runExPANdS(\"", snv_file, "\",\"", cn_file, "\",maxScore=2.5, max_PM=6, min_CellFreq=0.1, precision=NA,plotF=2,snvF=\"out.expands\",maxN=8000,region=NA)\n"]) )    
                    info_fh.write( "".join(["setwd(\"",wkdir,"\")\n"]))

def make_sankey_input_file( patient_files ):
    for patient in patient_files:
        #sps = sub populations
        primary_d = dict()
        refractory_d = dict()
        primary_snvd = dict()
        refractory_snvd = dict()
        #wkdir = "/projects/trans_scratch/validations/workspace/szong/IF-AML/genome/expands/test/" 
        wkdir = os.getcwd() + "/"
        sankey_header = ["primary_sp", "mutation", "refractory_sp", "value" ]
        if (len(patient_files[patient]) == 3):
             all_mutations = []
             primary_sps = "".join([wkdir, patient, "_diagnosis.expands.sps"])
             primary_snv= "".join([wkdir, patient, "_diagnosis.expands.snv"])
             refractory_sps = "".join([wkdir, patient, "_unknow_mix.expands.sps" ])
             refractory_snv = "".join([wkdir, patient, "_unknow_mix.expands.snv" ])
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
                             #primary_sp = "_".join(["relapse",str(primary_sp)])
                     except KeyError:
                         primary_sp = "primary_NA"
                         #primary_sp = "relapse_NA"
                     try:
                         refractory_sp = refractory_d[mutation][0]
                         if ("NA" not in refractory_sp):
                             refractory_sp = "{:.3f}".format(float(refractory_sp)/max_refractory_sp)
                             #refractory_sp = "_".join(["refractory",str(refractory_sp)] )
                             refractory_sp = "_".join(["postmortem",str(refractory_sp)] )
                     except KeyError:
                         #refractory_sp = "refractory_NA"
                         refractory_sp = "postmortem_NA"
                     try:
                         geneX = primary_snvd[mutation][0]
                     except KeyError:
                         geneX = refractory_snvd[mutation][0]
                     mutation = "_".join([geneX, mutation])
                     sankey_width = 0.075
                     if ("MODIFIER" not in mutation):
                         sankey_width = 0.2
                     writer.writerow([primary_sp, mutation, refractory_sp, sankey_width] )



if __name__ == '__main__':
    __main__()


