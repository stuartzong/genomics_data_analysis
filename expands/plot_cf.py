
import os
import subprocess
import re
import sys
import glob
#import argparse
import csv
from collections import defaultdict
from pprint import pprint
"""
This script does the following:

"""

def make_variant_patient_dict( snv_summary_file, patient_CRstatus, cf_change_file, outfile ):
    variant_patients = dict()
    with open (snv_summary_file, 'r') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        for line in records:
            patient_status = line['PatientID']
            gene = line['Gene']
            chr = line['Chromosome']
            pos = line['position']
            ref = line['referenceBase']
            alt = line['AlternativeBase']
            #variant = "_".join([ gene, chr, pos, ref, alt ])
            variant = "_".join([ gene, chr, pos ])
            try:
                CR_status = patient_CRstatus[ patient_status ]
            except KeyError:
                CR_status = "CR_status_unknown"
            value = "_".join([patient_status, CR_status  ])
            try:
                #print patient, status
                variant_patients[ variant ].append( value)
            except KeyError:
                    variant_patients[ variant ] = [ value ]



    with open ( outfile,  'wb' ) as writer:
        tree = lambda: defaultdict(tree)
        geneSnvPatientDict = tree()
        geneDict = dict()
        with open (cf_change_file, 'r') as fh:
            records = csv.DictReader(fh,  delimiter='\t')
            headers = records.fieldnames
            final_headers = ["gene", "numberKnownCRstatusPatientsCGIgeneLevel", "numberTypeOfCellularFrequencyChangeGeneLevel", "numberSensitveOrresistantSNVs"] + headers + ["numberOfPatientsInAllCGIgenomes", "patientsWithKnownCRstatusInCGIgenomes"]
            writer.write( "\t".join(final_headers))
            writer.write("\n")

            for line in records:
                sl = line['mutation'].split("_")
                gene = sl[0]
                type = line['TypeOfMutation']
                variant = "_".join(sl[i] for i in [0, 2, 3])
                #print "oooo",variant
                try: 
                    patients_all = variant_patients[ variant ]
                    patients_CRstatus_known = [i for i in patients_all if "CR_status_unknown" not in i ]
                    num_all_pats = len(patients_all) 
                    num_CRstatus_known_pats = len(patients_CRstatus_known)
                    spc_patients = ";".join(patients_CRstatus_known ) 
                except KeyError:
                    spc_patients = "variant_not_in_SPC_patient"
                    patients_CRstatus_known = [] 
                    num_all_pats = 0
                try:
                      geneDict[gene].extend( patients_CRstatus_known )
                except:
                      geneDict[gene] = patients_CRstatus_known

                content = [line[i] for i in headers ]
                #print content
                content.extend([ str(num_all_pats), spc_patients ])
                details = ":".join( content )
                try:
                    geneSnvPatientDict[gene][type][variant].append(details)
                except:
                    geneSnvPatientDict[gene][type][variant] = [details]
        print "length, ", len(geneSnvPatientDict)
        for gene in geneSnvPatientDict:
            numPatientGeneLevel = len(list(set(geneDict[gene])))
            numTypes = len (geneSnvPatientDict[gene].keys())
            for type in geneSnvPatientDict[gene]:
                numSNVs = int(len(geneSnvPatientDict[gene][type]))
                for snv in geneSnvPatientDict[gene][type]:
                    for content in geneSnvPatientDict[gene][type][snv]:
                        details = content.split(":")
                        #details = list(set(geneSnvPatientDict[gene][type][snv]))[0].split(":")
                        writer.write("\t".join([gene,  str(numPatientGeneLevel), str(numTypes), str(numSNVs),
                                                details[0],  details[1],  details[2],  details[3],  details[4], str(details[5]), str(details[6])]))
                                                
                        #writer.write("\t".join(content))
                        writer.write("\n")
                        print "ppppp", gene, type, snv, str(numPatientGeneLevel), str(numTypes)
        #pprint( geneSnvPatientDict )
        return variant_patients
    


def make_patient_CRstatus_dict( clinical_file ):
    #patient -> status -> CR status
    patient_CRstatus = dict()
    with open (clinical_file, 'r') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        for line in records:
            patient = line['TARGET USI'].split("-")[2] 
            status1 = "Primary"
            status2 = line['First Event']
            CRstatus1 = line['CR status at end of course 1']
            CRstatus2 = line['CR status at end of course 2']
            key1 = "_".join([ patient, status1 ]) 
            key2 = "_".join([ patient, status2 ]) 
            patient_CRstatus[ key1 ] = CRstatus1
            patient_CRstatus[ key2 ] = CRstatus2
        #pprint (patient_CRstatus)
        return patient_CRstatus



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

            try:
                #print patient, status
                patient_status_files[patient][status] = value
            except KeyError:
                #print "key error!"
                if patient not in patient_status_files:
                    patient_status_files[patient] = {}
                if status not in patient_status_files[patient]:
                    patient_status_files[patient][status] = value
        #print patient_status_files
        return patient_status_files



def count_SNVs( variant_summary ):
   #outfile = ".".join([ variant_summary, "filtered" ])
   outfile = "SNV_counts.txt"
   headers = ["patient", "numberOfAllMutationsPrimary", "numberOfNonslienceMutationsPrimary", "numberOfAllMutationsRefractory", "numberOfNonslienceMutationsRefractory"]
   d = dict()
   with open ( outfile,  'wb' ) as writer:
       #writer = csv.writer( fh, delimiter='\t' )

       print "Reading variant summary file is: %s.\n" % variant_summary
       with open ( variant_summary, 'r' ) as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            #headers = records.fieldnames
            writer.write( "\t".join(headers))
            writer.write("\n")
            for line in records:
                gene = line["Gene"]
                chr = line["Chromosome"]
                pos = line["position"]
                ref = line["referenceBase"]
                alt = line["AlternativeBase"]
                pat_status = line["PatientID"].split("_")
                pat = pat_status[0]
                status = pat_status[1]
                variant = "_".join([gene, chr, pos, ref, alt])
                #details_d = {key: line[key] for key in headers}
                try:
                    d[pat][status].append(variant)
                except KeyError:
                    if pat not in d:
                        d[pat] = {}
                    if status not in d[pat]:
                        d[pat][status] = [variant]

            for pat in d:
                count = []
                for status in d[pat]:
                    #print pat, status
                    types = ["HIGH", "MODERATE"]
                    variants = list(set( d[pat][status] ))
                    num_all_variants = len( variants )
                    #nonsilence_variants = [i for i in variants if any(type in variants for type in types )]
                    nonsilence_variants = [i for i in variants if any( type in i for type in types) ]
                    num_nonsilence_variants = len( nonsilence_variants )
                    count.extend([ num_all_variants, num_nonsilence_variants])
                tmp = "\t".join([str(i) for i in count])     
                count = [i for i in count] 
                #columns = [  d[variant][pat][header] for header in headers]
                writer.write( "\t".join([pat, tmp]))
                writer.write("\n")
       

def make_sps_dict( patient_files ):
    #sps_d = cancer cell subpopulations dictionary
    #cf = cellular frequency
    header = ["primary", "gene", "impact", "refractory"]
    for patient in patient_files:
        wkdir = os.getcwd()
        wkdir = wkdir + "/"
        if (len(patient_files[patient]) == 3):
             all_mutations = []
             sankey = "".join([wkdir, patient, ".sankey"])
             out_file = "".join([wkdir, patient, "_cf_plot.csv"])
             with open ( out_file,  'wb' ) as fh3:
                writer = csv.writer( fh3, delimiter='\t' )
                writer.writerow( header )
                with open ( sankey, 'r' ) as fh1:
                     recd1 = csv.DictReader(fh1,  delimiter='\t')
                     for line in recd1:
                         primary_sp = line["primary_sp"].split("_")[1]
                         if ("NA" in primary_sp):
                             primary_sp = 0
                         refractory_sp = line["refractory_sp"].split("_")[1]
                         if ("NA" in refractory_sp):
                             refractory_sp = 0
                         sl = line["mutation"].split("_")
                         if (len(sl) == 4):
                             gene = sl[0]
                             impact = sl[1]
                             writer.writerow([str(primary_sp), gene, impact, str(refractory_sp)])
    
            
def __main__():

    print "\nGenerating patient_files dictionary!"
    input_files = "bam_vcf_cnv_path.txt"
    patient_files = make_patient_vcf_file_dict(input_files)
    pprint(patient_files)

    make_sps_dict(patient_files)

if __name__ == '__main__':
    __main__()


