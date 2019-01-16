import os
import subprocess
import re
import sys
import glob
#import argparse
import csv
from collections import defaultdict


def parse_maf(CGI_variant_details_file, targetted_genes):
    gene_variant_dict = dict()
    with open(CGI_variant_details_file, 'rU') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        out_file = "CGI_variant_details_file_new.txt"
        with open(out_file, 'wb') as writer:
            for line in records:
                try:
                    #print line
                    var_type  = line['VariantType']
                    if (var_type == ""):
                        var_type = "NA"
                    mut_status = line["Mutation_Status"]
                    if (mut_status == ""):
                        mut_status = "NA"
                    chr = line['Chromosome'].strip()
                    start = line['Start_position'].strip()
                    end = line['End_position'].strip()
                    gene = line['Hugo_Symbol']
                    if (gene == ""):
                        gene = "NA"
                    var_classification = line['Variant_Classification']
                    if (var_classification == ""):
                        var_classification = "NA"
                    dbsnp = line['dbSNP_RS'] 
                    if (dbsnp == ""):
                        dbsnp = "NA"
                    tum_alt = line['Tumor_ReadCount_Alt']	
                    tum_ref = line['Tumor_ReadCount_Ref']	
                    tum_total = line['Tumor_ReadCount_Total']	
                    nor_alt = line['Normal_ReadCount_Alt']	
                    nor_ref = line['Normal_ReadCount_Ref']	
                    nor_total = line['Normal_ReadCount_Total']	
                    cosm_id = line['Cosmic']	
                    if (cosm_id == ""):
                        cosm_id = "NA"
                    cosm_gene = line['Cosmic_Gene']	
                    if (cosm_gene == ""):
                        cosm_gene = "NA"
                    ref = line['Reference_Allele']	
                    tum_allele1 = line['TumorSeq_Allele1']	
                    tum_allele2 = line['TumorSeq_Allele2']	
                    if (tum_allele2 == ""):
                        tum_allele2 = "N"
                    if (tum_allele1 == ref):
                        alt = tum_allele2
                    else:
                        alt = tum_allele1
                    nor_allele1 = line['Match_Norm_Seq_Allele1']	
                    nor_allele2 = line['Match_Norm_Seq_Allele2']	
                    if (nor_allele2 == ""):
                        nor_allele2 = "N"
                    tum_barcode = line['Tumor_Sample_Barcode']
                    variant = ":".join([chr, start, ref, alt])
                    #patient = "-".join(tum_barcode.split("-")[0:3])
                    patient = tum_barcode.split("-")[2]
                    value = ":".join([patient, nor_total, nor_ref, nor_alt, tum_total, tum_ref, tum_alt])
                    #if (mut_status == "Somatic") and (var_type == "SNP") and (gene in targetted_genes):
                    if ((var_type == "SNP") and (gene in targetted_genes)):
                        #dict can have duplicate patients when a snp is seen in both primary and relapse
                        try:
                            if value not in gene_variant_dict[gene][variant]:
                                gene_variant_dict[gene][variant].append(value)
                        except KeyError:
                            if gene not in  gene_variant_dict:
                                gene_variant_dict[gene] = {}
                            if variant not in gene_variant_dict[gene]:
                                gene_variant_dict[gene][variant] = [value]
                        writer.write("\t".join([chr, start, gene, ref, alt]))
                        writer.write("\n") 
                except:#there are apparent mistakes in the CGI maf files
                    print "Invalid entry in original CGI maf file! This line is skipped!\n", 
                    continue
        print "CGI gene_variant_dict", gene_variant_dict 
        return gene_variant_dict   
    
def list_files_by_extension(extension, vcfFileList):
    filelist = [ f for f in os.listdir(".") if f.endswith(tuple(extension))]
    writer = open(vcfFileList, 'w')
    for f in filelist:
        #print "writing vcf files into vcfFileList: %s" % f
        writer.write(f)
        writer.write("\n")


def analyze_pileup_fisher_file(variant_gene_dict, fisher_files): 
    fish_gene_variant_dict = {}
    with open(fisher_files, 'r') as fh:
        for f in fh:
            fish_file = f.strip()
            with open(fish_file,  'r') as handle:
                patient = "_".join(re.split("[-|.]", fish_file)[2:4])
                handle.next()
                for line in handle:
                    #print line
                    sl = line.strip().split("\t")
                    chr, start, ref, alt = sl[1].strip(), sl[2].strip(), sl[3].strip(), sl[4].strip()
                    if (alt == "U"):
                        alt = sl[9].strip()
                    nor_total, nor_ref = sl[5].strip(), sl[6].strip()
                    nor_alt, nor_af = sl[7].strip(), sl[8].strip()
                    tum_total, tum_ref = sl[10].strip(), sl[11].strip(), 
                    tum_alt, tum_af = sl[12].strip(), sl[13].strip()
                    adj_tum_total, adj_tum_ref = sl[14].strip(), sl[15].strip(), 
                    adj_tum_alt, adj_tum_af = sl[16].strip(), sl[17].strip()
                    tum_con = sl[18].strip()
                    try:
                        p_value = float(sl[19].strip())
                    except:
                        p_value = 1
                    #filter for homozygous reference in normal OR heterzygous in normal and lost of normal copy in tumur
                    if ((int(nor_total) > 10 and int(tum_total) > 10 and int(tum_alt) > 3 and int(nor_alt) < 2) or (
                        int(nor_total) > 10 and int(tum_total) > 10 and float(adj_tum_af) >0.8 and float(nor_af) <0.6
                        and float(nor_af) > 0.1)) and p_value < 0.05 : 
   
                        variant = ":".join([chr, start, ref, alt])
                        #if pileup was run on more positions, subset to postions we are trying to verify
                        try: 
                            gene = variant_gene_dict[variant]
                            value = ":".join([patient, nor_total, nor_ref, nor_alt, nor_af,
                                              tum_total, tum_ref, tum_alt, tum_af,
                                              adj_tum_total, adj_tum_ref, adj_tum_alt, adj_tum_af,
                                              tum_con, str(p_value)])
                            try:
                                if value not in fish_gene_variant_dict[gene][variant]:
                                    fish_gene_variant_dict[gene][variant].append(value)
                            except KeyError:
                                if gene not in  fish_gene_variant_dict:
                                    fish_gene_variant_dict[gene] = {}
                                if variant not in fish_gene_variant_dict[gene]:
                                    fish_gene_variant_dict[gene][variant] = [value]
                        except KeyError:
                            print "%s is not one of the variants to be verified!" % variant
                            continue 
    print "fish_gene_variant_dict", fish_gene_variant_dict
    return fish_gene_variant_dict

def make_dis_val_patient_list(infile):
    dis_val_patients = []
    with open(infile) as dis_handle:
        for line in dis_handle:
            #patient_name = "-".join(line.strip().split("-")[0:4])
            patient_name = line.strip().split("-")[2]
            dis_val_patients.append(patient_name)
    #print "dis_val_patients are: %s\t%s\n" % (len(dis_val_patients), dis_val_patients)
    no_dup_patients = list(set(dis_val_patients))
    print "number of patient in both discovery and validation is: %s\n%s\n" % (len(no_dup_patients), no_dup_patients) 
    return dis_val_patients

def make_targetted_gene_list(targetted_gene_file):
    targetted_genes = []
    with open(targetted_gene_file) as tar_handle:
        for line in tar_handle:
            gene_name = line.strip()
            targetted_genes.append(gene_name)
    print "targetted_genes are: %s\n%s\n" % (len(targetted_genes), targetted_genes)
    return targetted_genes


def make_spc_som_variant_dict(spc_som_genes_summary, targetted_genes):
    spc_variant_dict = dict()
    with open (spc_som_genes_summary,  'r') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        for line in records:
            gene = line["Gene"]
            patient = line["PatientID"]
            chr = line["Chromosome"]
            pos = line["position"]
            ref = line["referenceBase"]
            alt = line["AlternativeBase"]
            variant = "_".join([chr, pos, ref, alt])
            if (gene in targetted_genes):

                try:
                    spc_variant_dict[gene][variant].append(patient)
                except KeyError:
                    print "keyerror\n"
                    if gene not in  spc_variant_dict:
                        spc_variant_dict[gene] = {}
                    if variant not in spc_variant_dict[gene]:
                        spc_variant_dict[gene][variant] = [patient]
        print spc_variant_dict
        return spc_variant_dict

def make_variant_list(val_gene_dict, spc_variant_dict):
        variants_list = []
        for gene in val_gene_dict:
            for variant in val_gene_dict[gene]:

                #for variant in variants:
                    gene_variant = "_".join([gene, variant])
                    variants_list.append(gene_variant)
                    print "variant is %s\n" % gene_variant
        for gene in spc_variant_dict:
            for variant in spc_variant_dict[gene]:
                #for variant in variants:
                    gene_variant = "_".join([gene, variant])
                    variants_list.append(gene_variant)
                    print "variant is %s\n" % gene_variant
        variants_list = list(set(variants_list))
        print variants_list
        return variants_list

def compare_2_nested_dicts(val_gene_dict, spc_som_gene_dict, variants_list, out_file):
    with open(out_file, 'wb') as writer:
        writer.write("\t".join(["gene", "variant", "numberOfPatientInDiscovery", "patientInDiscovery", "numberOfPatientInValidation", "patientInValidation", "discoveryDetails", "validationDetails"]))
        writer.write("\n")
        for variant in variants_list:
            print variant
            sp = variant.split("_")
            gene = sp[0]
            variant = "_".join(sp[1:])
            try:
                pat_discovery = val_gene_dict[gene][variant]
                num_discovery = len(pat_discovery)
            except KeyError:
                num_discovery = 0
                pat_discovery = "N"
           
            try:
                pat_spc = spc_som_gene_dict[gene][variant]
                num_spc = len(pat_spc)
                
                pat_shared = set(pat_discovery) & set(pat_spc)
                #print "pat_shared:%s%s%s\n" % (pat_spc, dis_val_patients, pat_shared)
                num_shared = len(pat_shared)
                #num_addi = num_spc - num_shared
                #pat_addi= set(pat_spc) - set(pat_shared)
                if num_shared == 0:
                    pat_shared = "N"
                #if num_addi ==0:
                #    pat_addi = "N"
                
            except KeyError:
                num_spc = 0
                pat_spc = "N"
                num_shared = 0
                pat_shared = "N"
                '''
                num_addi = 0
                pat_addi = "N"
                

                try:
                    pat_cgi = cgi_gen_som_gene_dict[gene]
                    num_cgi = len(pat_cgi)
                except KeyError:
                    pat_cgi = "N"
                    num_cgi = 0
                '''
            print gene, variant, num_discovery, pat_discovery, num_spc, pat_spc, num_shared, pat_shared
            
            writer.write("\t".join([gene, variant, str(num_discovery), ", ".join( pat_discovery), str(num_spc), ", ".join( pat_spc), str(num_shared), ", ".join( pat_shared) ]))
            writer.write("\n")
            

def compare_nested_dicts(CGI_gene_variant_dict, fish_gene_variant_dict, out_file):
    #make full gene:variant list
    with open(out_file, 'wb') as writer:
        writer.write("\t".join(["gene", "variant", "numberOfPatientInDiscovery", "patientInDiscovery", "numberOfPatientInValidation", "patientInValidation", "discoveryDetails", "validationDetails"]))
        writer.write("\n")
        all_gene_variants = []
        for gene in CGI_gene_variant_dict:
            for variant in CGI_gene_variant_dict[gene]:
                gene_var = "_".join([gene, variant])
                all_gene_variants.append(gene_var)
        for gene in fish_gene_variant_dict:
            for variant in fish_gene_variant_dict[gene]:
                gene_var = "_".join([gene, variant])
                all_gene_variants.append(gene_var)
        all_gene_variants = list(set(all_gene_variants))
        print "all_gene_variants", all_gene_variants
        
        for gv in all_gene_variants:
            sl = gv.split("_")
            gene, variant = sl[0], sl[1] 
            try:
                CGI_details = CGI_gene_variant_dict[gene][variant]
            except KeyError:
                CGI_details = ["None"]
            try:
                fish_details = fish_gene_variant_dict[gene][variant]
            except KeyError:
                fish_details = ["None"]
            CGI_patients = list(set([i.split(":")[0] for i in CGI_details]))
            num_CGI_patients = str(len(CGI_patients))
            fish_patients = [i.split(":")[0] for i in fish_details]
            num_fish_patients = len(list(set([i.split("_")[0] for i in fish_patients])))
            if (fish_details == ["None"]):
                num_fish_patients-=1
            print gene, variant, len(CGI_patients), CGI_patients, len(fish_patients), fish_patients, CGI_details, fish_details
        
            writer.write("\t".join([gene, variant, num_CGI_patients, ";".join(CGI_patients), str(num_fish_patients),
                         ";".join(fish_patients), ";".join(CGI_details), ";".join(fish_details)]))   
            writer.write("\n") 
def __main__():
    wkdir = "/projects/trans_scratch/validations/workspace/szong/AML_capture/"
    targetted_gene_file = "".join([wkdir, "meta/targetted_genes.txt"])
    CGI_variant_details_file = "/projects/trans_scratch/validations/workspace/szong/AML_capture/patients_in_both/fullbarcode_matched/sqhigh_and_verified_variants_in_genes_somatic_non_silent_details_with_filters_AML_May_14_2014.tsv" 
    #CGI_variant_details_file = "86_patient_sqhigh_and_verified_variants_in_genes_somatic_non_silent_details_with_filters_AML_May_14_2014.tsv"

    #targetted gene list, 396 genes
    targetted_genes = make_targetted_gene_list(targetted_gene_file)

    #make cgi verified gene>variant dict exluding untargeted genes
    CGI_gene_variant_dict = parse_maf(CGI_variant_details_file, targetted_genes)      
    '''
    # make a mapping dict gene <> variants
    variants_gene_mapping_dict ={}
    for gene in CGI_gene_variant_dict:
        for variant in CGI_gene_variant_dict[gene]:
            variants_gene_mapping_dict[variant] = gene
    print variants_gene_mapping_dict

    #list all fisher files
    extension = ['combined.adjusted.fisher']
    fisher_files = "fisher.files"
    list_files_by_extension(extension, fisher_files)
     



    #analyze pileup fisher files
    fish_gene_variant_dict = analyze_pileup_fisher_file(variants_gene_mapping_dict, fisher_files)
    

    #compare 2 dict
    out_file = "validated_variants_comparison.txt"
    compare_nested_dicts(CGI_gene_variant_dict, fish_gene_variant_dict, out_file)
    '''
 
if __name__ == '__main__':
    __main__()

