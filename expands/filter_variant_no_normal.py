
import os
import subprocess
import re
import sys
import glob
#import argparse
import csv
from collections import defaultdict

"""
This script does the following:

"""
def make_dbsnp142_flagged_list(flagged_file):
   flagged = []
   with open ( flagged_file, 'r' ) as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            for line in records:
                rs = line["name"]
                flagged.append(rs)
   return flagged


def make_gene_list(gene_file):
   genes = []
   with open ( gene_file, 'r' ) as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            for line in records:
                gene = line["name"]
                genes.append(gene)
   return genes

def make_gene_rank(rank_file):
   gene_ranks = {}
   with open ( rank_file, 'r' ) as handle:
        for line in handle:
            sl = line.strip().split("\t")
            gene = sl[1]
            rank = sl[0]
            if (gene not in gene_ranks):
                gene_ranks[gene] = rank
   return gene_ranks


def filter_variants(somatic_summary):
   filtered_summary = ".".join([somatic_summary, "filtered" ])
   """ d = gene_variant_patients dictionary """
   d = dict()
   diagnosis_variants = []
   relapse_variants = []
   pm_variants = []
   diagnosis_file = 'diagnosis_variants.txt'
   relapse_file = 'relapse_variants.txt'
   pm_file = 'postmortem_variants.txt'
   with open (filtered_summary,  'wb') as fh:
       writer = csv.writer( fh, delimiter='\t' )
       print "The variant summary file is: %s.\n" % somatic_summary
       with open (somatic_summary, 'r') as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            headers = records.fieldnames
            writer.writerow(headers)
            for line in records:
                gene = line['gene']
                chr = line["chromosome"]
                pos = line["position"]
                ref = line["referenceBase"]
                alt = line["alternativeBase"]
                DNA_cov = int(line["tumourDNASequencingCoverage"])
                DNA_refC = int(line["tumourDNARefBaseCount"])
                DNA_altC = int(line["tumourDNAAltBaseCount"])
                DNA_af = float(line["tumourDNAAlleleFrequency"])
                RNA_cov = int(line["tumourRNASequencingCoverage"])
                RNA_refC = int(line["tumourRNARefBaseCount"])
                RNA_altC = int(line["tumourRNAAltBaseCount"])
                RNA_af = float(line["tumourRNAAlleleFrequency"])
                RNA_altref_total = RNA_refC + RNA_altC
                content = [line[i] for i in headers]
                pass_filter = False
                if ((DNA_cov > 10 and DNA_altC > 3 and DNA_af >= 0.1) or
                   (RNA_cov > 10 and RNA_altC > 3 and RNA_af >= 0.1 and RNA_altref_total > 0.8*RNA_cov)): 
                    writer.writerow(content)
                    patient = line['patientID']
                    variant =":".join([gene, "_".join([chr, pos, ref, alt])])
                    if ('diagnosis' in  patient):
                        print "diagnosis", patient 
                        diagnosis_variants.append(variant)
                    elif ('relapse' in  patient):
                        print "relapse", patient 
                        relapse_variants.append(variant)
                    elif ('mix' in  patient):
                        print "pm", patient 
                        pm_variants.append(variant)
                    #print "True var''iant\n"
                    #print "\t".join(content)
                else:
                    #print "fake variant\n"
                    print "\t".join(content)
       write_list(diagnosis_variants, diagnosis_file)
       write_list(relapse_variants, relapse_file)
       write_list(pm_variants, pm_file)

def write_list(lst, file):                     
    with open (file, 'wb') as writer:
        variants = list(set(lst))
        for variant in variants:
            sl = re.split('[:_]', variant)
            writer.write("\t".join(sl))
            writer.write('\n')

def count_patients( filtered_summary, genvar_pat_dict):
   final_summary = filtered_summary + ".final"
   header = ["Rank", "Gene", "numPatientGeneLevel", "numVariantsForThisGene", 
              "numPatientGeneLevelFiltered", "numVariantsForThisGeneFiltered", 
              "Chromosome", "position", "referenceBase", "AlternativeBase", 
              "numPatientVariantLevel", "numPatientVariantLevelFiltered",
                "PatientID", "dbSNPid", "clinic", "globalMinorAlleleFrequency", "COSMIC64ID",
                "ImpactDetails", "normalSequencingCoverage", "normalRefBaseCount",
                "normalAltBaseCount", "normalAlleleFrequency",
                "tumourSequencingCoverage", "tumourRefBaseCount",
                "tumourAltBaseCount", "tumourAlleleFrequency",
                "adjustedTumourSequencingCoverage","adjustedTumourRefBaseCount",
                "adjustedTumourAltBaseCount", "adjustedTumourAlleleFrequency",
                "fisherExactPvalue", "tumourContent", "cosmicAnnotationV71",
                "inStrelka"]
   with open ( final_summary,  'wb' ) as fh:
       writer = csv.writer( fh, delimiter='\t' )
       writer.writerow( header )
       with open ( filtered_summary, 'r' ) as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            for line in records:
                rank = line["Rank"]
                gene = line["Gene"]
                clinic = line["clinic"]
                numpat_gene = line["numPatientGeneLevel"]
                num_var = line["numVariantsForThisGene"]
                chr = line["Chromosome"]
                pos = line["position"]
                ref = line["referenceBase"]
                alt = line["AlternativeBase"]
                numpat_var = line["numPatientVariantLevel"]
                pat = line["PatientID"]
                dbsnp = line["dbSNPid"]
                gmaf = line["globalMinorAlleleFrequency"]
                cos_id = line["COSMIC64ID"]
                impact = line["ImpactDetails"]
                t_cov = line["tumourSequencingCoverage"]
                t_refC = line["tumourRefBaseCount"]
                t_altC = line["tumourAltBaseCount"]
                t_af = line["tumourAlleleFrequency"]
                tc = line["tumourContent"]
                cos_anno = line["cosmicAnnotationV71"]
                shorted = line["shortlisted"]

                variant = "_".join([chr, pos, ref, alt])
                gene_vars = list(set([i for i in genvar_pat_dict[gene].keys()]))
                num_var_filtered = len(gene_vars)
                tmp_gene_pats = []
                for variant1 in gene_vars:
                    tmp_gene_pats = tmp_gene_pats + genvar_pat_dict[gene][variant1]
                gene_pats = [i.split("_")[0] for i in tmp_gene_pats]               
                numpat_gene_filtered = len(list(set(gene_pats))) 
                var_pats = [i.split("_")[0] for i in genvar_pat_dict[gene][variant]]
                numpat_var_filtered = len(list(set(var_pats)))
                writer.writerow([ rank, gene, numpat_gene, num_var, numpat_gene_filtered, num_var_filtered, 
                                  chr, pos, ref, alt, numpat_var, numpat_var_filtered, 
                                  pat, dbsnp, clinic, gmaf, cos_id, impact,
                                          t_cov, t_refC, t_altC, t_af,
                                          tc, cos_anno, shorted])

def make_list(infile):
   variants = []
   with open (infile, 'r' ) as fh:
        for line in fh:
            sl = line.rstrip().split('\t')
            variant = "_".join(sl)
            variants.append(variant)
   return variants

def filter_somatic(variant_summary, somatic_variants):
   filtered_summary = ".".join([ variant_summary, "filtered" ])
   with open (filtered_summary,  'wb') as fh:
       writer = csv.writer(fh, delimiter='\t')
       print "The variant summary file is: %s.\n" % variant_summary
       with open (variant_summary, 'r') as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            headers = records.fieldnames
            for line in records:
                content = [line[i] for i in headers] 
                chr = line["chromosome"]
                pos = line["position"]
                ref = line['referenceBase']
                alt = line['alternativeBase']
                variant = "_".join([chr, pos, ref, alt])
                if (variant in somatic_variants):
                    writer.writerow(content) 


def filter_variants_333(somatic_summary):
   filtered_summary = ".".join([somatic_summary, "filtered" ])
   """ d = gene_variant_patients dictionary """
   d = dict()
   right_variants = []
   left_variants = []
   right_file = '333_right_variants.txt'
   left_file = '333_left_variants.txt'
   with open (filtered_summary,  'wb') as fh:
       writer = csv.writer( fh, delimiter='\t' )
       print "The variant summary file is: %s.\n" % somatic_summary
       with open (somatic_summary, 'r') as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            headers = records.fieldnames
            writer.writerow(headers)
            for line in records:
                gene = line['gene']
                chr = line["chromosome"]
                pos = line["position"]
                ref = line["referenceBase"]
                alt = line["alternativeBase"]
                DNA_cov = int(line["tumourDNASequencingCoverage"])
                DNA_refC = int(line["tumourDNARefBaseCount"])
                DNA_altC = int(line["tumourDNAAltBaseCount"])
                DNA_af = float(line["tumourDNAAlleleFrequency"])
                RNA_cov = int(line["tumourRNASequencingCoverage"])
                RNA_refC = int(line["tumourRNARefBaseCount"])
                RNA_altC = int(line["tumourRNAAltBaseCount"])
                RNA_af = float(line["tumourRNAAlleleFrequency"])
                RNA_altref_total = RNA_refC + RNA_altC
                content = [line[i] for i in headers]
                pass_filter = False
                if ((DNA_cov > 10 and DNA_altC > 3 and DNA_af >= 0.1) or
                   (RNA_cov > 10 and RNA_altC > 3 and RNA_af >= 0.1 and RNA_altref_total > 0.8*RNA_cov)): 
                    writer.writerow(content)
                    patient = line['patientID']
                    variant =":".join([gene, "_".join([chr, pos, ref, alt])])
                    if ('right' in  patient):
                        print "right", patient
                        right_variants.append(variant)
                    elif ('left' in  patient):
                        print "left", patient
                        left_variants.append(variant)
                else:
                    #print "fake variant\n"
                    print "\t".join(content)
       write_list(right_variants, right_file)
       write_list(left_variants, left_file)

                         


def filter_variants_hvsn(somatic_summary):
   filtered_summary = ".".join([somatic_summary, "filtered" ])
   """ d = gene_variant_patients dictionary """
   d = dict()
   hypo_variants = []
   normoxic_variants = []
   hypo_file = '333_hypo_variants.txt'
   normoxic_file = '333_normoxic_variants.txt'
   with open (filtered_summary,  'wb') as fh:
       writer = csv.writer( fh, delimiter='\t' )
       print "The variant summary file is: %s.\n" % somatic_summary
       with open (somatic_summary, 'r') as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            headers = records.fieldnames
            writer.writerow(headers)
            for line in records:
                gene = line['gene']
                chr = line["chromosome"]
                pos = line["position"]
                ref = line["referenceBase"]
                alt = line["alternativeBase"]
                DNA_cov = int(line["tumourDNASequencingCoverage"])
                DNA_refC = int(line["tumourDNARefBaseCount"])
                DNA_altC = int(line["tumourDNAAltBaseCount"])
                DNA_af = float(line["tumourDNAAlleleFrequency"])
                RNA_cov = int(line["tumourRNASequencingCoverage"])
                RNA_refC = int(line["tumourRNARefBaseCount"])
                RNA_altC = int(line["tumourRNAAltBaseCount"])
                RNA_af = float(line["tumourRNAAlleleFrequency"])
                RNA_altref_total = RNA_refC + RNA_altC
                content = [line[i] for i in headers]
                pass_filter = False
                if ((DNA_cov > 10 and DNA_altC > 3 and DNA_af >= 0.1) or
                   (RNA_cov > 10 and RNA_altC > 3 and RNA_af >= 0.1 and RNA_altref_total > 0.8*RNA_cov)): 
                    writer.writerow(content)
                    patient = line['patientID']
                    variant =":".join([gene, "_".join([chr, pos, ref, alt])])
                    if ('hypo' in  patient):
                        print "hypo", patient
                        hypo_variants.append(variant)
                    elif ('normoxic' in  patient):
                        print "normoxic", patient
                        normoxic_variants.append(variant)
                else:
                    #print "fake variant\n"
                    print "\t".join(content)
       write_list(hypo_variants, hypo_file)
       write_list(normoxic_variants, normoxic_file)



def __main__():
   #variant_summary = "high_moderate_SNV_summary_no_normal.txt"
   #potenital_somatics = "/projects/trans_scratch/validations/workspace/szong/David_Kaplan/variants/somatic/all_potential_somatic_snps.txt"
   #somatic_variants = make_list(potenital_somatics)
   #print somatic_variants
   #filter_somatic(variant_summary, somatic_variants)
   #filtered_somatic = "high_moderate_SNV_summary_no_normal.txt.128153295.somatic"
   #filter_variants(filtered_somatic)

   # filter 333 normoxic right vs left
   
   #rl_333_normoxic_somatic = "high_moderate_SNV_summary_no_normal.txt.333.RL.normoxic.somatic"
   #filter_variants_333(rl_333_normoxic_somatic)
   #filter 333 right hypo vs normoxic
   #right_333_hvsn_somatic = "high_moderate_SNV_summary_no_normal.txt.333.right.somatic"
   #filter_variants_hvsn(right_333_hvsn_somatic)

   left_333_hvsn_somatic = "high_moderate_SNV_summary_no_normal.txt.333.left.somatic"
   filter_variants_hvsn(left_333_hvsn_somatic)







if __name__ == '__main__':
    __main__()


