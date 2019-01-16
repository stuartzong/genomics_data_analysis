
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


def filter_variants( variant_summary, flagged_snp, genes, gene_ranks ):
   filtered_summary = ".".join([ variant_summary, "filtered" ])
   """ d = gene_variant_patients dictionary """
   d = dict()
   header = ["Rank", "Gene", "numPatientGeneLevel",
                "numVariantsForThisGene", "Chromosome", "position", "referenceBase",
                "AlternativeBase", "numPatientVariantLevel",
                "PatientID", "dbSNPid", "clinic", "globalMinorAlleleFrequency", "COSMIC64ID",
                "ImpactDetails",
                "tumourSequencingCoverage", "tumourRefBaseCount", 
                "tumourAltBaseCount", "tumourAlleleFrequency",
                "tumourContent", "cosmicAnnotationV71", 
                "shortlisted"]
   with open ( filtered_summary,  'wb' ) as fh:
       writer = csv.writer( fh, delimiter='\t' )
       writer.writerow( header )

       print "The variant summary file is: %s.\n" % variant_summary
       with open ( variant_summary, 'r' ) as handle:
            records = csv.DictReader(handle,  delimiter='\t')
            for line in records:
                gene = line["Gene"]
                if (gene in gene_ranks):
                    rank = gene_ranks[gene]
                else:
                    rank = "unknown"
                numpat_gene = int(line["numPatientGeneLevel"])
                num_var = line["numVariantsForThisGene"]
                chr = line["Chromosome"]
                pos = line["position"]
                ref = line["referenceBase"]
                alt = line["AlternativeBase"]
                numpat_var = int(line["numPatientVariantLevel"])
                pat = line["PatientID"]
                dbsnp = line["dbSNPid"]
                try:
                    gmaf = float(line["globalMinorAlleleFrequency"])
                except:#include gmaf unknown calls
                    gmaf = 0
                cos_id = line["COSMIC64ID"]
                impact = line["ImpactDetails"]
                t_cov = int(line["tumourSequencingCoverage"])
                t_refC = line["tumourRefBaseCount"]
                t_altC = int(line["tumourAltBaseCount"])
                t_af = float(line["tumourAlleleFrequency"])
                tc = line["tumourContent"]
                cos_anno = line["cosmicAnnotationV71"]
                clinic = "clinic_unknown"
                shorted = "not_shortlisted"
                pass_filter = False
                if (t_cov > 10 and t_altC > 3 and
                    t_af >= 0.05 and gmaf < 0.01): 
                    if (dbsnp in flagged_snp):
                        clinic = "clinic_associated"
                        pass_filter = True
                    # <10% total patients if gene in cosmic
                    elif ( numpat_var < 15 and (("not_in_COSMIC71" not in cos_anno) or
                           ("not_in_COSMIC" not in cos_id)) ):
                        pass_filter = True
                    #<10% patient if gene is in the gene short list
                    elif ((gene in genes)  and
                           numpat_var > 0  and numpat_var < 15):
                        shorted = "shortlisted"
                        pass_filter = True
                    #otherwise, occurrence < 3% of total patients
                    elif ( numpat_var > 0  and numpat_var < 5):
                        pass_filter = True
                    
                    t_af = line["tumourAlleleFrequency"]
                    gmaf = line["globalMinorAlleleFrequency"]     
                    t_altC = line["tumourAltBaseCount"]
                    if(pass_filter):
                        writer.writerow([ rank, gene, numpat_gene, num_var, chr, pos, ref, alt,
                                          numpat_var, pat, dbsnp, clinic, gmaf, cos_id, impact,
                                          t_cov, t_refC, t_altC, t_af,
                                          tc, cos_anno, shorted])
                        variant = "_".join([chr, pos, ref, alt])
                        try:
                            d[gene][variant].append( pat )
                        except KeyError:
                            #print "key error!"
                            if gene not in d:
                                d[gene] = {}
                            if variant not in d[gene]:
                                d[gene][variant] = [ pat ]   
       print d
       return d


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

def filter_variants(variant_summary, somatic_variants):
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
                         
def __main__():
   variant_summary = "high_moderate_SNV_summary_no_normal.txt"
   potenital_somatics = "/projects/trans_scratch/validations/workspace/szong/David_Kaplan/variants/somatic/all_potential_somatic_snps.txt"
   somatic_variants = make_list(potenital_somatics)
   print somatic_variants
   filter_variants(variant_summary, somatic_variants)
   sys.exit()
   interested_genes = "knownAMLmutatedGenes.txt"
   genes = make_gene_list(interested_genes)
   rank_file = "gene_rank.txt"
   gene_ranks = make_gene_rank(rank_file)

   flagged_file = "snp142Flagged.txt"
   flagged_snp = make_dbsnp142_flagged_list(flagged_file)
   print flagged_snp
   genvar_pat_dict = filter_variants( variant_summary, flagged_snp, genes, gene_ranks)
   filtered_summary = ".".join([ variant_summary, "filtered" ])
   count_patients(filtered_summary, genvar_pat_dict)
if __name__ == '__main__':
    __main__()


