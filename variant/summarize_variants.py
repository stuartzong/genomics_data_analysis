#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice

# test comments
"""
testing script: expands pipeline (v1.0) for production
To run this script, type the following command:
python  variant_summarization_clonality_analysis.py -dt wgs or spc > summary_log_file.txt
normal: matched normal samples must be specified as "normal", case sensitive
"""
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


def take(n, iterable):
    """ Return first n items of the iterable as a list! """
    return list(islice(iterable, n))

def is_other_readable(file_path):
    m = os.stat(file_path).st_mode
    #otherExec  = bool(m & 0001)
    #otherWrite = bool(m & 0002)
    #other_read_status  = bool(m & 0004)
    group_read_status = bool(m & stat.S_IRGRP)
    return group_read_status
    #return other_read_status


def exist(file_path):
    exist_status = os.path.exists(file_path)
    return exist_status
 
def check_files(files):
    missing_files = []
    if (len(files) == 0):
       print "No files to check!"
    for file in files:
        short_name = file.split("/")[-1]
        exist_status = False
        read_status = False
        exist_status = exist(file)
        if (exist_status):
            read_status = is_other_readable(file)
        if (not(exist_status and read_status)):
            missing_files.append(file)
        else:
            print "%s: OK" % short_name 
    if missing_files:
        print "ERROR: The following files are either missing or no read permission! Check if STRELKA has completed!"
        for file in missing_files:
            print "Missing %s\n" % file
            sys.exit()

def make_pileup_script(patient, status, data_type, pos_file, bam_file):
    """ Generate mpileup scripts! """
    reference_genome= "/home/pubseq/genomes/Homo_sapiens/hg19a/bwa_ind/genome/GRCh37-lite.fa"
    samtools_path = "/home/rcorbett/aligners/samtools/samtools-0.1.17/samtools"
    python_path = "/gsc/software/linux-x86_64/python-2.7.2/bin/python"
    python_script = "/home/szong/python/PileupParser_CountSnpBases.py"
    bash_file = ".".join([patient, status, data_type, "pileup.sh"])
    pileup_outfile = ".".join([patient, status, data_type, "pileup"])
    base_count_file  = ".".join([patient, status, data_type, "pileup.AFcounts"])
    with open(bash_file,  'wb') as writer:
         writer.write("#! /bin/bash\n")
         writer.write("#$ -S /bin/bash\n")
         writer.write("#$ -N samtools\n")
         writer.write("#$ -q transabyss.q\n")
         writer.write("#$ -l mem_token=2G,mem_free=2G,h_vmem=2G\n")
         writer.write("#$ -V\n")
         writer.write(" ".join([samtools_path, "mpileup -f", reference_genome,
                                "-l", pos_file, bam_file, ">", pileup_outfile, "\n"]))
         #writer.write(" ".join(["rm", pos_file, "\n"]))
         writer.write("if [ $? -eq 0 ]\n")
         writer.write("    then\n")
         writer.write("    touch %s_%s_%s_pileup.complete\n" % (patient, status, data_type)) 
         writer.write(" ".join(["echo", "\"mpileup job finished successfully at:\" `date`", "\n"])) 
         writer.write("else\n")
         writer.write("    echo \"ERROR: samtools mpileup did not finish correctly!\" >>\"summary_log_file.txt\"\n")
         writer.write("    echo \"Failed patient is: %s_%s_%s \" >> \"summary_log_file.txt\"\n" % (patient, status, data_type)) 
         writer.write("    exit\n")
         writer.write("fi\n")
         writer.write(" ".join([python_path, python_script, "-f -i ",
                                pileup_outfile, "-o", base_count_file, "\n"]))
         writer.write(" ".join(["echo", "\"mpileup and post-processing job finished successfully at:\" `date`", "\n"])) 
    return bash_file


def delete_files_by_extension(extension):
    files = [ f for f in os.listdir(".") if f.endswith(tuple(extension))]
    for f in files:
        os.remove(f)

def delete_files(files):
    for f in files:
        os.remove(f)

def run_fisher_exact_test(Rscript_path,r_script):
    p = subprocess.Popen('%s/Rscript %s' % (Rscript_path,r_script),  shell=True, stdout=subprocess.PIPE)
    output,  err = p.communicate()
    #print output

def count_ref_alt_bases(snp_list, base_count_file, af_file):
   """ Generate allele frequency dict to append af info to summary file """ 
   snps = dict()
   for snp in snp_list:
       sl = snp.rstrip().split(':')
       chr = sl[0] 
       pos = sl[1]
       ref = sl[2]
       alt = sl[3]
       key = ":".join([chr, pos, ref])
       if key not in snps:
           snps[key] = [alt]
       else:
           snps[key].append(alt)

   # parse mpileup output to make af dict
   afs = dict() 
   print "Parsing mpileup file: %s!" % base_count_file
   with open (base_count_file) as fh:
       for line in fh:
           sl = line.strip().split('\t')
           pos_ref = sl[0].split(":")
           chr = pos_ref[0]
           pos = pos_ref[1]
           ref = pos_ref[2]
           snp = sl[0]
           cov = sl[1] 
           a, c, g, t, n = sl[2], sl[3], sl[4], sl[5], sl[6]
           d = {"A":a, "C":c, "G":g, "T":t, "N":n}
           ref_count = d[ref]
           d.pop(ref, None)
           alts = snps[snp]
           for alt in alts:
               alt_count = d[alt]
               adj_cov = int(alt_count)+int(ref_count)
               if (adj_cov != 0):
                   af = str("{:.2f}".format(float(alt_count)/adj_cov))
               else:
                   af = str("{:.2f}".format(float(alt_count)/int(cov)))
               value = ":".join([alt, cov, ref_count, alt_count, af])

               try:
                   afs[snp].append(value)
               except KeyError:
                   afs[snp] = [value]

   # print afs
   with open(af_file,  'wb') as opf:
       writer = csv.writer(opf,  delimiter='\t')
       for snp in afs:
           sp_snp = snp.split(':')
           for af in afs[snp]:
               sp_af = af.split(':')
               content = sp_snp + sp_af
               writer.writerow(content)
   return afs


def parse_vcf(vcf_file, caller, impacts, gene_pats, gene_variants, snv_pos, snv_list, patient_status):
    print "%s, Variant caller is: %s!\n " %  (vcf_file.split("/")[-1], caller)
    
    with open(vcf_file,  'r') as fh:
        ## Following line is useful for testing the script
        #fh = [next(fh) for x in xrange(350)]
        for line in fh:
            if (not line.startswith("#")):
                sl = line.strip().split("\t")
                Info = sl[7]
                if (not Info.startswith("INDEL")) and ("EFF=" in Info) and any(x in Info for x in impacts):
                    splitInfo = Info.split(";")
                    # get gmaf value      
                    try:
                       gmaf = [j for j in splitInfo if
                           j.startswith("GMAF=")][0].split("=")[1]
                    except:
                       gmaf = "gmaf_unknown"
                    transcripts = [i for i in splitInfo if
                        i.startswith("EFF=")][0].split("=")[1].split(",")
                    rs_cos_id = sl[2]
                    chr, start, ref = sl[0], sl[1], sl[3]
                    alts = sl[4].split(",")
                    snv_pos.append("\t".join([chr, start]))
                    for alt in alts:
                        snv = ":".join([chr, start, ref, alt])
                        snv_list.append(snv)
                        snp_id_tmp = rs_cos_id.split(";")[0]
                        if snp_id_tmp.startswith("rs"):
                           snp_id = snp_id_tmp
                        else:
                           snp_id = "novel_snp"
                        try:
                           cosmic_id = rs_cos_id.split(";")[1]
                        except:
                           if snp_id_tmp.startswith("COS"):
                               cosmic_id = snp_id_tmp
                           else:
                               cosmic_id = "not_in_cosmic64"
                        selected_transcripts = [k for k in transcripts if any(x in k for x in impacts)]
                        # pick the first transcript to get the gene name
                        trans_details = selected_transcripts[0]
                        gene = trans_details.split("|")[5]
                        if (gene == ""):
                            gene = "INTERGENIC"
                        # make gene -> patient dict to count how many patients have variant in this gene
                        try:
                           gene_pats[gene].append(patient_status)
                        except:
                           gene_pats[gene] = [patient_status]
                        details = ":".join([snp_id, cosmic_id, gmaf, trans_details, caller])
                        # put all relevant info into a triple nested dictionary: gene -> snv -> patient -> details
                        try:
                           gene_variants[gene][snv][patient_status].append(details)
                        except:
                           gene_variants[gene][snv][patient_status] = [details]
    return [gene_pats, gene_variants, snv_pos, snv_list]
    

def summarize_snvs(patient_files, sum_header_wn_tmp, hm_sum_wn_tmp, hm_impacts):
    """ 
    Generate triple nested dictionary for each
    Combine all the SNV positions from each pileline and patient status for pileup 
    Positions = mpileup + strelka + primary + relaspe
    Get snp list for all patients
    Summarize all SNVs with high or moderate or low or modifier impact
    Annoated as if a SNV in mpileup or/and strelka
    pool all snv positions from various vcf files: single, paired, mutseq, strelka
    if a variant exists in one tissue, time point, check all samples from the same patient
    """
    tree = lambda: defaultdict(tree)
    hm_variants_dict = tree()
    #low_variants_dict = tree()
    hm_gene_dict = dict()
    #low_gene_dict = dict()
    patient_snv = dict()
    patient_callers = dict()
    snv_pos_files = []
    #hm_impacts = ["HIGH", "MODERATE", "LOW"]
    for patient in patient_files:
        print "Start parsing all vcf files related to patient: %s\n" % patient
        snv_pos_file = ".".join([patient, "vcf.snp.pos"]) 
        snv_pos_files.append(snv_pos_file)
        snv_pos = []        
        snv_list = []
        for status in patient_files[patient]:
            if ("normal" not in status):
                print "Parsing vcf files related to patient and status: %s, %s\n" % (patient, status)
                strelka_snv_vcf = patient_files[patient][status]['strelka_snv_vcf']
                DNA_single_vcf = patient_files[patient][status]['DNA_single_vcf']
                paired_mpileup_vcf = patient_files[patient][status]['paired_mpileup_vcf']
                mutseq_snv_vcf = patient_files[patient][status]['mutseq_snv_vcf']
                RNA_single_vcf = patient_files[patient][status]['RNA_single_vcf']
                patient_status = "_".join([patient, status])
    
                if (strelka_snv_vcf != "NA"):
                    caller = "strelka"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_vcf(strelka_snv_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, snv_pos, snv_list, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]
                    snv_pos = two_dicts[2]
                    snv_list = two_dicts[3]
                if (mutseq_snv_vcf != "NA"):
                    caller = "mutseq"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_vcf(mutseq_snv_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, snv_pos, snv_list, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]
                    snv_pos = two_dicts[2]
                    snv_list = two_dicts[3]
    
                if (DNA_single_vcf != "NA"):
                    caller = "DNA_single_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_vcf(DNA_single_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, snv_pos, snv_list, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]
                    snv_pos = two_dicts[2]
                    snv_list = two_dicts[3]
    
                if (paired_mpileup_vcf != "NA"):
                    caller = "paired_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_vcf(paired_mpileup_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, snv_pos, snv_list, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]
                    snv_pos = two_dicts[2]
                    snv_list = two_dicts[3]

                if (RNA_single_vcf != "NA"):
                    caller = "RNA_single_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_vcf(RNA_single_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, snv_pos, snv_list, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]
                    snv_pos = two_dicts[2]
                    snv_list = two_dicts[3]

        snv_pos = list(set(snv_pos))
        snv_list = list(set(snv_list))
        patient_snv[patient] = snv_list
        with open(snv_pos_file, 'w') as pos_writer:
            for snv_coord in sorted(snv_pos):
                pos_writer.write(snv_coord)
                pos_writer.write("\n")
 
    # combine varians in both mpileup and strelka vcfs
    # print patient,"\n", hm_variants_dict
    #print "Example patient snv dict are: \n"
    #n_items = take(2, patient_snv.iteritems())
    #print n_items
    write_summary(hm_variants_dict, hm_gene_dict, hm_sum_wn_tmp, sum_header_wn_tmp, patient_callers, patient_files)
    return [patient_snv, snv_pos_files]
    
def write_summary(variants_dict, gene_dict, summary_file, sum_header_wn_tmp, patient_callers, patient_files):
    #write info into summary file
    writer = open(summary_file, 'w')
    writer.write("\t".join(sum_header_wn_tmp))
    writer.write("\n")
    print "Start writing variants_dict into a summary file: %s\n" % summary_file
    for gene in variants_dict:
        short_patient = [i.split("_")[0] for i in gene_dict[gene]]
        num_patients_gene_level = len(list(set(short_patient)))
        numSNVs = len(variants_dict[gene].keys())
        for snv in variants_dict[gene]:
            short_num_patients_SNV_level = [i.split("_")[0] for i in variants_dict[gene][snv]]
            num_patients_SNV_level = len(list(set(short_num_patients_SNV_level)))
            snv_details = []
            combination_in_snv = [i for i in variants_dict[gene][snv]]
            #print "short_patients in patient_files: %s\n" % patient_file
            #print "short_num_patients_SNV_level are: %s\n" % short_num_patients_SNV_level
            # all patient status combinations
            patients_snv_level = list(set(short_num_patients_SNV_level))
            all_combinations = []
            for pat in patients_snv_level: 
                    for sta in patient_files[pat]:
                        if ("normal" not in sta):
                            pat_sta = "_".join([pat, sta])
                            all_combinations.append(pat_sta)
            #print "aaaaaa", patients_snv_level
            combination_not_called = list(set(all_combinations)-set(combination_in_snv))
            #print "all_combinations is %s, %s, %s, %s\n" % (gene, all_combinations, combination_in_snv, combination_not_called)

            for patient in variants_dict[gene][snv]:
                #print "patient and status_num_patients_SNV_level is %s, %s\n" % (patient, status_num_patients_SNV_level)
                # add relapse in if only primary in variants_dict
                # snv_details: 10:64927823:C:G
                snv_details = snv.split(":")
                chr = snv_details[0]
                start = snv_details[1]
                ref = snv_details[2]
                alt = snv_details[3]
                #anno_details: rs8473:COSM1128106:0.4748:NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aaa/Gaa|K2857E|2896|MKI67|protein_coding|CODING|ENST00000368653|13|1):pileup
                anno_details_all = variants_dict[gene][snv][patient]
                anno_details = list(set(anno_details_all))[0].split(":")
                callers = [i.split(":")[-1] for i in anno_details_all]
                rs_id = anno_details[0]
                gmaf = anno_details[2]
                cosmic_id = anno_details[1]
                snp_eff = anno_details[3]
                #print "writing gene, snv, patient,anno_details_all: %s, %s, %s, %s\n" % (gene, snv, patient,anno_details_all)
                # see if variant called by both mpileup and strelka
                in_paired_mpileup = "not_in_paired_mpileup"
                in_mutseq = "not_in_mutseq"
                in_strelka = "not_in_strelka"
                in_DNA_single_mpileup = "not_in_DNA_single_mpileup"
                in_RNA_single_mpileup = "not_in_RNA_single_mpileup"
                if ("strelka" in patient_callers[patient]):
                    if ("strelka" in callers):
                        in_strelka = "in_strelka"
                else:
                    in_strelka = "strelka_not_run"

                if ("paired_mpileup" in patient_callers[patient]):
                    if ("paired_mpileup" in callers):
                        in_paired_mpileup = "in_paired_mpileup"
                else:
                    in_paired_mpileup = "paired_mpileup_not_run"

                if ("DNA_single_mpileup" in patient_callers[patient]):
                    if ("DNA_single_mpileup" in  callers):
                        in_DNA_single_mpileup = "in_DNA_single_mpileup"
                else:
                    in_DNA_single_mpileup = "DNA_single_mpileup_not_run"

                if ("RNA_single_mpileup" in patient_callers[patient]):
                    if ("RNA_single_mpileup" in  callers):
                        in_RNA_single_mpileup = "in_RNA_single_mpileup"
                else:
                    in_RNA_single_mpileup = "RNA_single_mpileup_not_run"

                if ("mutseq" in patient_callers[patient]):
                    if ("mutseq" in callers):
                        in_mutseq = "in_mutseq"
                else:
                    in_mutseq = "mutseq_not_run"

                content = "\t".join([gene, str(num_patients_gene_level), str(numSNVs),
                                        chr, start, ref, alt, str(num_patients_SNV_level), 
                                        patient, rs_id, gmaf, cosmic_id, snp_eff, in_DNA_single_mpileup,
                                        in_paired_mpileup, in_mutseq, in_strelka, in_RNA_single_mpileup])
                writer.write(content)
                writer.write("\n")
            for combination in combination_not_called:    
                patient = combination
                # see if variant called by both mpileup and strelka
                in_paired_mpileup = "not_in_paired_mpileup"
                in_mutseq = "not_in_mutseq"
                in_strelka = "not_in_strelka"
                in_DNA_single_mpileup = "not_in_DNA_single_mpileup"
                in_RNA_single_mpileup = "not_in_RNA_single_mpileup"
                if ("strelka" not in patient_callers[patient]):
                    in_strelka = "strelka_not_run"
 
                if ("paired_mpileup" not in patient_callers[patient]):
                    in_paired_mpileup = "paired_mpileup_not_run"
 
                if ("DNA_single_mpileup" not in patient_callers[patient]):
                    in_DNA_single_mpileup = "DNA_single_mpileup_not_run"
 
                if ("RNA_single_mpileup" not in patient_callers[patient]):
                    in_RNA_single_mpileup = "RNA_single_mpileup_not_run"

                if ("mutseq" not in patient_callers[patient]):
                    in_mutseq = "mutseq_not_run"

                content = "\t".join([gene, str(num_patients_gene_level), str(numSNVs),
                                        chr, start, ref, alt, str(num_patients_SNV_level), 
                                        combination, rs_id, gmaf, cosmic_id, snp_eff, in_DNA_single_mpileup,
                                        in_paired_mpileup, in_mutseq, in_strelka, in_RNA_single_mpileup])

                writer.write(content)
                writer.write("\n")



def calculate_af (cov, refC, altC):
    # ref + alt <= cov because low qaulity bases filtered
    adj_cov = int(refC) + int(altC)
    if (adj_cov != 0):
        af = "{:.2f}".format(float(altC)/adj_cov)
    elif (int(cov) != 0):
        af = "{:.2f}".format(float(altC)/int(cov))
    else:
        af = 0.00
    return af


def summarize_indels(patient_files, sum_header_wn_tmp, hm_sum_wn_tmp, hm_impacts):
    tree = lambda: defaultdict(tree)
    hm_variants_dict = tree()                                                                                         
    #low_variants_dict = tree()
    modifier_variants_dict = tree()
    hm_gene_dict = dict()
    #low_gene_dict = dict()
    modifier_gene_dict = dict()
    hm_indel_dict = dict()
    #low_indel_dict = dict()
    patient_callers = dict()    
    modifier_indel_dict = dict()
    #hm_impacts = ["HIGH", "MODERATE"]
    #low_impacts = ["LOW"]

    for patient in patient_files:
        print "Start parsing all vcf files related to patient: %s\n" % patient
        for status in patient_files[patient]:
            if ("normal" not in status):
                print "Parsing all vcf files related to patient and status: %s, %s\n" % (patient, status)
                strelka_indel_vcf = patient_files[patient][status]['strelka_indel_vcf']
                DNA_single_vcf = patient_files[patient][status]['DNA_single_vcf']
                RNA_single_vcf = patient_files[patient][status]['RNA_single_vcf']
                paired_mpileup_vcf = patient_files[patient][status]['paired_mpileup_vcf']
                #mutseq does not have indel result
                #mutseq_snv_vcf = patient_files[patient][status]['mutseq_snv_vcf']
                patient_status = "_".join([patient, status])

                if (strelka_indel_vcf != "NA"):
                    caller = "strelka"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_strelka_indel(strelka_indel_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]

                if (DNA_single_vcf != "NA"):
                    caller = "DNA_single_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_mpileup_indel(DNA_single_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]

                if (RNA_single_vcf != "NA"):
                    caller = "RNA_single_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_mpileup_indel(RNA_single_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]


                if (paired_mpileup_vcf != "NA"):
                    caller = "paired_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_mpileup_indel(paired_mpileup_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]

    write_indel_summary(hm_variants_dict, hm_gene_dict, hm_sum_wn_tmp, sum_header_wn_tmp, patient_callers, patient_files)


def parse_mpileup_indel(vcf_file, caller, impacts, gene_pats, gene_variants, patient_status):                                                          
    
    print "%s, Variant caller is: %s!\n " %  (vcf_file.split("/")[-1], caller)
    with open(vcf_file,  'r') as fh:
        ## Following line is useful for testing the script
        # pileup header 153 lines
        #fh = [next(fh) for x in xrange(200)]
        for line in fh:
            if (not line.startswith("#")):
                sl = line.strip().split("\t")
                Info = sl[7]
                if (Info.startswith("INDEL")) and ("EFF=" in Info) and any(x in Info for x in impacts):
                    #print Info, impacts
                    splitInfo = Info.split(";")
                    cov =  [i for i in splitInfo if i.startswith("DP=")][0].split("=")[1]
                    refC_altC =  [i for i in splitInfo if i.startswith("DP4=")][0].split("=")[1].split(",")
                    tmp_refC = refC_altC[0:2]
                    tmp_altC = refC_altC[2:]
                    refC = sum(int(i) for i in tmp_refC)
                    altC = sum(int(i) for i in tmp_altC)
                    af = calculate_af (cov, refC, altC)     
                    mpileup_cov = [cov, str(refC), str(altC), str(af)]
                    strelka_cov = ["na", "na", "na", "na", "na", "na", "na", "na"]
                    cov_info = "\t".join(mpileup_cov + strelka_cov)
                    try:
                        gmaf = [j for j in splitInfo if
                                j.startswith("GMAF=")][0].split("=")[1]
                    except:
                       gmaf = "gmaf_unknown"
                    transcripts = [i for i in splitInfo if
                        i.startswith("EFF=")][0].split("=")[1].split(",")
                    rs_cos_id = sl[2]
                    chr, start, ref = sl[0], sl[1], sl[3]
                    alts = sl[4].split(",")
                    for alt in alts:
                        indel = ":".join([chr, start, ref, alt])
                        snp_id_tmp = rs_cos_id.split(";")[0]
                        if snp_id_tmp.startswith("rs"):
                           snp_id = snp_id_tmp
                        else:
                           snp_id = "novel_snp"

                        try:
                           cosmic_id = rs_cos_id.split(";")[1]
                        except:
                           if snp_id_tmp.startswith("COS"):
                               cosmic_id = snp_id_tmp
                           else:
                               cosmic_id = "not_in_cosmic64"

                        selected_transcripts = [k for k in transcripts if any(x in k for x in impacts)]
                        # pick the first transcript to get the gene name
                        trans_details = selected_transcripts[0]
                        gene = trans_details.split("|")[5]
                        if (gene == ""):
                            gene = "INTERGENIC"
                        # make gene -> patient dict to count how many patients have variant in this gene
                        try:
                           gene_pats[gene].append(patient_status)
                        except:
                           gene_pats[gene] = [patient_status]
                        details = ":".join([snp_id, cosmic_id, gmaf, trans_details, cov_info, caller])
                        # put all relevant info into a triple nested dictionary: gene -> snv -> patient -> details
                        try:
                           gene_variants[gene][indel][patient_status].append(details)
                        except:
                           gene_variants[gene][indel][patient_status] = [details]
    return [gene_pats, gene_variants]


    

def parse_strelka_indel(vcf_file, caller, impacts, gene_pats, gene_variants, patient_status):                                                          
    print "%s, Variant caller is: %s!\n " %  (vcf_file.split("/")[-1], caller)
    with open(vcf_file,  'r') as fh:
        ## Following line is useful for testing the script
        #strelka header lines 324
        #fh = [next(fh) for x in xrange(300)]
        for line in fh:
            if (not line.startswith("#")):
                sl = line.strip().split("\t")
                Info = sl[7]
                if ("EFF=" in Info) and any(x in Info for x in impacts):
                            # use tier1 counts, which is with mapping quality >40
                            normal_info = sl[9].split(":")
                            tumor_info = sl[10].split(":")
                            splitInfo = Info.split(";")
                            nor_cov =  normal_info[0]
                            nor_refC = normal_info[2].split(",")[0]
                            nor_indelC = normal_info[3].split(",")[0]
                            tum_cov =  tumor_info[0]
                            tum_refC = tumor_info[2].split(",")[0]
                            tum_indelC = tumor_info[3].split(",")[0]
                            #nor_adj_cov = int(nor_refC) + int(nor_indelC)
                            #tum_adj_cov = int(tum_refC) + int(tum_indelC)
                            nor_af = calculate_af (nor_cov, nor_refC, nor_indelC)
                            tum_af = calculate_af (tum_cov, tum_refC, tum_indelC)
                            mpileup_cov = ["na", "na", "na", "na"]
                            strelka_cov = [nor_cov, nor_refC, nor_indelC, str(nor_af), tum_cov, tum_refC, tum_indelC, str(tum_af)]
                            cov_info = "\t".join(mpileup_cov + strelka_cov)
    
                            try:
                                gmaf = [j for j in splitInfo if
                                        j.startswith("GMAF=")][0].split("=")[1]
                            except:
                               gmaf = "gmaf_unknown"
                            transcripts = [i for i in splitInfo if
                                i.startswith("EFF=")][0].split("=")[1].split(",")
                            rs_cos_id = sl[2]
                            chr, start, ref = sl[0], sl[1], sl[3]
                            alts = sl[4].split(",")
                            for alt in alts:
                                indel = ":".join([chr, start, ref, alt])
                                snp_id_tmp = rs_cos_id.split(";")[0]
                                if snp_id_tmp.startswith("rs"):
                                   snp_id = snp_id_tmp
                                else:
                                   snp_id = "novel_snp"
    
                                try:
                                   cosmic_id = rs_cos_id.split(";")[1]
                                except:
                                   if snp_id_tmp.startswith("COS"):
                                       cosmic_id = snp_id_tmp
                                   else:
                                       cosmic_id = "not_in_cosmic64"

                                selected_transcripts = [k for k in transcripts if any(x in k for x in impacts)]
                                trans_details = selected_transcripts[0]
                                gene = trans_details.split("|")[5]
                                if (gene == ""):
                                    gene = "INTERGENIC"
                                try:
                                   gene_pats[gene].append(patient_status)
                                except:
                                   gene_pats[gene] = [patient_status]
                                details = ":".join([snp_id, cosmic_id, gmaf, trans_details, cov_info, caller])
                                try:
                                   gene_variants[gene][indel][patient_status].append(details)
                                except:
                                   gene_variants[gene][indel][patient_status] = [details]
    return [gene_pats, gene_variants]





def write_indel_summary(variants_dict, gene_dict, summary_file, sum_header_wn_tmp, patient_callers, patient_files):
    #write info into summary file
    writer = open(summary_file, 'w')
    writer.write("\t".join(sum_header_wn_tmp))
    writer.write("\n")
    print "Start writing variants_dict into a summary file: %s\n" % summary_file
    for gene in variants_dict:
        short_patient = [i.split("_")[0] for i in gene_dict[gene]]
        num_patients_gene_level = len(list(set(short_patient)))
        numSNVs = len(variants_dict[gene].keys())
        for snv in variants_dict[gene]:
            short_num_patients_SNV_level = [i.split("_")[0] for i in variants_dict[gene][snv]]
            num_patients_SNV_level = len(list(set(short_num_patients_SNV_level)))
            snv_details = []
            for patient in variants_dict[gene][snv]:
                # snv_details: 10:64927823:C:G
                snv_details = snv.split(":")
                chr = snv_details[0]
                start = snv_details[1]
                ref = snv_details[2]
                alt = snv_details[3]
                #anno_details: rs8473:COSM1128106:0.4748:NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|Aaa/Gaa|K2857E|2896|MKI67|protein_coding|CODING|ENST00000368653|13|1):pileup
                anno_details_all = variants_dict[gene][snv][patient]
                #anno_details = list(set(anno_details_all))[0].split(":")
                if ("strelka" in anno_details_all): 
                    anno_details = [i for i in anno_details_all if "strelka" in i][0].split(":")
                else:
                    anno_details = anno_details_all[0].split(":") 
                callers = [i.split(":")[-1] for i in anno_details_all]
                rs_id = anno_details[0]
                gmaf = anno_details[2]
                cosmic_id = anno_details[1]
                snp_eff = anno_details[3]
                cov_info = anno_details[4]
                # see if variant called by both mpileup and strelka
                in_paired_mpileup = "not_in_paired_mpileup"
                in_strelka = "not_in_strelka"
                in_DNA_single_mpileup = "not_in_DNA_single_mpileup"
                in_RNA_single_mpileup = "not_in_RNA_single_mpileup"
                if ("strelka" in patient_callers[patient]):
                    if ("strelka" in callers):
                        in_strelka = "in_strelka"
                else:
                    in_strelka = "strelka_not_run"

                if ("paired_mpileup" in patient_callers[patient]):
                    if ("paired_mpileup" in callers):
                        in_paired_mpileup = "in_paired_mpileup"
                else:
                    in_paired_mpileup = "paired_mpileup_not_run"

                if ("DNA_single_mpileup" in patient_callers[patient]):
                    if ("DNA_single_mpileup" in  callers):
                        in_DNA_single_mpileup = "in_DNA_single_mpileup"
                else:
                    in_DNA_single_mpileup = "DNA_single_mpileup_not_run"

                if ("RNA_single_mpileup" in patient_callers[patient]):
                    if ("RNA_single_mpileup" in  callers):
                        in_RNA_single_mpileup = "in_RNA_single_mpileup"
                else:
                    in_RNA_single_mpileup = "RNA_single_mpileup_not_run"

                content = "\t".join([gene, str(num_patients_gene_level), str(numSNVs),
                                        chr, start, ref, alt, str(num_patients_SNV_level), 
                                        patient, rs_id, gmaf, cosmic_id, snp_eff, in_DNA_single_mpileup,
                                        in_paired_mpileup, in_strelka, in_RNA_single_mpileup, cov_info])
                writer.write(content)
                writer.write("\n")
  

def adjust_allele_frequency(type, com_af_files,patient_files):
    tc_key = "_".join([type, "tc"])
    for file in com_af_files:
        file = file.strip()
        sp_f = re.split('[.]',file)
        patient =sp_f[0]
        status = sp_f[1]
        pat_status = "_".join([patient, status])
        af_aj_file = file + ".adjusted"
        writer = csv.writer(open(af_aj_file, 'wb'),delimiter='\t')
        with open (file) as handle:
            for line in handle:
                sl = line.rstrip().split('\t')
                n_totC = int(sl[4])
                n_refC = int(sl[5])
                n_altC = int(sl[6])
                t_totC = int(sl[9])
                t_refC = int(sl[10])
                t_altC = int(sl[11])
                n_af = float(sl[7])
                t_af = float(sl[12])
                # skip adjustment for second low frequency allele 
                tlow_totC = t_refC + t_altC
                tc_tmp = patient_files[patient][status][tc_key]
                if ( n_totC > 0 and t_totC > 0 and (tlow_totC > 0.7*t_totC) and ("unknow" not in tc_tmp) ):
                    nlz_factor = float(float(t_totC)/n_totC)
                    tc = float(int(tc_tmp)/100.0)
                    n_nlz_refC = n_refC*nlz_factor
                    n_nlz_altC = n_altC*nlz_factor
                    t_aj_ref = float(t_refC-n_nlz_refC*(1-tc))/tc
                    t_aj_alt = float(t_altC-n_nlz_altC*(1-tc))/tc
                    if (t_aj_ref >= 0 and t_aj_alt >= 0):
                        t_aj_tot = t_aj_ref + t_aj_alt
                        if t_aj_tot != 0:
                            t_aj_af = t_aj_alt/t_aj_tot
                            if t_aj_af >= 1:
                                t_aj_af =1
                            elif t_aj_af<0:
                                t_aj_af = 0
                        else:
                            t_aj_tot, t_aj_ref = t_totC, t_refC
                            t_aj_alt, t_aj_af = t_altC, t_af
                    elif (t_aj_ref < 0 and t_aj_alt >= 0):
                        t_aj_alt = t_aj_alt - t_aj_ref
                        t_aj_ref = 0
                        t_aj_tot = t_aj_ref + t_aj_alt
                        if t_aj_tot != 0:
                            t_aj_af = t_aj_alt/t_aj_tot
                            if t_aj_af >= 1:
                                t_aj_af =1
                            elif t_aj_af < 0:
                                t_aj_af = 0
                        else:
                            t_aj_tot, t_aj_ref = t_totC, t_refC
                            t_aj_alt, t_aj_af = t_altC, t_af
                    elif (t_aj_ref >= 0 and t_aj_alt < 0):
                        t_aj_ref = t_aj_ref - t_aj_alt
                        t_aj_alt = 0
                        t_aj_tot = t_aj_ref + t_aj_alt
                        if t_aj_tot != 0:
                            t_aj_af = t_aj_alt/t_aj_tot
                            if t_aj_af >= 1:
                                t_aj_af = 1
                            elif t_aj_af < 0:
                                t_aj_af = 0
                        else:
                            t_aj_tot, t_aj_ref = t_totC, t_refC
                            t_aj_alt, t_aj_af = t_altC, t_af
                else:
                    t_aj_tot, t_aj_ref = t_totC, t_refC
                    t_aj_alt, t_aj_af = t_altC, t_af
                    try:
                        tc = float(int(tc_tmp)/100.0)
                    except:
                        tc = "tc_unknown"
                t_aj_tot = "{:.1f}".format( float(t_aj_tot) )
                t_aj_ref = "{:.1f}".format( float(t_aj_ref) )
                t_aj_alt = "{:.1f}".format( float(t_aj_alt) )
                t_aj_af = "{:.2f}".format( float(t_aj_af) )
                writer.writerow([sl[0],sl[1],sl[2],sl[3],sl[4],sl[5],sl[6],sl[7],sl[8],
                                t_totC, t_refC, t_altC, t_af, t_aj_tot,t_aj_ref,t_aj_alt,t_aj_af,tc])
        print "Finishing af adjustment for: %s" % (pat_status)

  

def make_patient_wn_af_dict(af_files):
   patient_af_dict = dict()
   for file in af_files:
       if (file != "NA"):
           patient = "_".join(re.split('[.]', file)[:2])
           with open(file) as f:
              next( f )
              for line in f:
                 sl = [i.strip() for i in line.rstrip().split('\t')]
                 chr, pos, ref, alt = sl[1], sl[2], sl[3], sl[4]
                 key = "_".join([patient, chr, pos, ref, alt])
                 n_cov, n_ref, n_alt, n_af = sl[5], sl[6], sl[7], sl[8]
                 t_cov, t_ref, t_alt, t_af = sl[10], sl[11], sl[12], sl[13]
                 n_af = "{:.2f}".format( float(n_af) )
                 t_af = "{:.2f}".format( float(t_af) )
                 ajt_cov, ajt_ref, ajt_alt, ajt_af = sl[14], sl[15], sl[16], sl[17]
                 ajt_cov = "{:.1f}".format( float(ajt_cov) )
                 ajt_ref = "{:.1f}".format( float(ajt_ref) )
                 ajt_alt = "{:.1f}".format( float(ajt_alt) )
                 ajt_af = "{:.2f}".format( float(ajt_af) )
                 p_value, tc = sl[19], sl[18]
               
                 value = ",".join([n_cov, n_ref, n_alt, n_af, 
                                   t_cov, t_ref, t_alt, t_af, 
                                   ajt_cov, ajt_ref, ajt_alt, ajt_af, 
                                   p_value])
                 try:
                   patient_af_dict[key].append(value)
                 except KeyError:
                   patient_af_dict[key] = [value]
   return patient_af_dict

def make_patient_non_af_dict(af_files):
   patient_af_dict = dict()
   for file in af_files:
       if (file != "NA"):
           patient = "_".join(re.split('[.]', file)[:2])
           with open(file) as f:
              for line in f:
                 sl = [i.strip() for i in line.rstrip().split('\t')]
                 chr, pos, ref, alt = sl[0], sl[1], sl[2], sl[3]
                 key = "_".join([patient, chr, pos, ref, alt])
                 t_cov, t_ref, t_alt, t_af = sl[4], sl[5], sl[6], sl[7]
                 value = ",".join([t_cov, t_ref, t_alt, t_af])
                 try:
                   patient_af_dict[key].append(value)
                 except KeyError:
                   patient_af_dict[key] = [value]
   return patient_af_dict



def qsub_scripts(scripts):
    """ qsub scripts """
    wkdir = os.getcwd()
    for script in scripts:
        #p = subprocess.Popen('qsub %s' % script,  shell=True, stdout=subprocess.PIPE)
        p = subprocess.Popen('ssh m0001 \"cd %s;  qsub %s\"' % (wkdir, script),  shell=True, stdout=subprocess.PIPE)
        output,  err = p.communicate()



def parse_pileup_output(patient_files, patient_snv_wn, patient_snv_non):
    af_files = dict()
    com_af_files = dict()
    af_dicts = dict()
    data_types = ["DNA", "RNA"]
    for type in data_types:
        af_files[type] = []
        com_af_files[type] = []
        af_dicts[type] = dict()
    for patient in patient_files:
        if ("normal" in patient_files[patient]):
            snp_list = patient_snv_wn[patient]
        else:
            snp_list = patient_snv_non[patient]

        # get base count  and af for all samples of a patient
        for status in patient_files[patient]:
            bam_types = dict()
            DNA_bam = patient_files[patient][status]['DNA_bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            bam_types["DNA"] = DNA_bam
            bam_types["RNA"] = RNA_bam
            for type in bam_types:
                if (bam_types[type] != "NA"):
                    #for type in data_types:
                    base_count_file  = ".".join([patient, status, type, "pileup.AFcounts"])
                    af_file = ".".join([patient, status, type, "pileup.AFcounts.af"])
                    af_files[type].append( af_file )
    
                    # variant bp count dict for primary, relapse, and normal
                    af_dict  = count_ref_alt_bases(snp_list, base_count_file, af_file)
                    # '5:131892979:G': ['G:88:1:1:0.01']
                    try:
                        #print patient, status
                        af_dicts[type][patient][status] = af_dict
                    except KeyError:
                        #print "key error!"
                        if patient not in af_dicts[type]:
                            af_dicts[type][patient] = {}
                        if status not in af_dicts[type][patient]:
                            af_dicts[type][patient][status] = af_dict
    
        nor_bams = dict()
        if ("normal" in patient_files[patient]):
            nor_bams["DNA"] = patient_files[patient]["normal"]['DNA_bam']
            nor_bams["RNA"] = patient_files[patient]["normal"]['RNA_bam']
        else:
            nor_bams["DNA"] = "NA"
            nor_bams["RNA"] = "NA"
        for type in nor_bams:    
            if (nor_bams[type] != "NA"):
                n_af_dict = af_dicts[type][patient]["normal"]
                for status in patient_files[patient]:
                    # if tumor bam exist
                    bam_key = "_".join([type, "bam"])
                    t_bam = patient_files[patient][status][bam_key]
                    if (not status == "normal") and (t_bam != "NA"):
                        #print "patient, status, bam_key, type are:", patient, status, bam_key, type
                        print "merging %s %s %s with its corresponding normal!" % (patient, status, type)
                        com_af_file = ".".join([patient, status, type, "af.combined"])
                        com_af_files[type].append(com_af_file)
                        t_af_dict = af_dicts[type][patient][status]
                        with open(com_af_file,  'wb') as opf:
                            writer = csv.writer(opf,  delimiter='\t')
                            #print "snp_list is: ", snp_list
                            for snp in snp_list:
                               sp_snp = snp.split(':')
                               chr, start, ref = sp_snp[0], sp_snp[1], sp_snp[2]
                               alt = sp_snp[3]
                               af_key = ":".join([chr, start, ref])
                               if ("no_matched_normal" not in n_af_dict):  
                                   try:
                                       normal_afs = n_af_dict[af_key]
                                       normal_af = [i for i in normal_afs if i[0] == alt][0].split(":")
                                       #print "normal_af is: %s, %s " % (type(normal_af), normal_af)
                                   except KeyError:
                                       normal_af = ["N", "0", "0", "0", "0"]
                                       #print "normal KEYERROR! %s \n" % normal_af
                               else:
                                   alt_tag = "_".join(["no", type, "normal"])
                                   normal_af = [alt_tag, "0", "0", "0", "0"]
                               n_totC, n_refC, n_altC, n_af = normal_af[1], normal_af[2], normal_af[3], normal_af[4]
    
                               try:
                                   tumour_afs = t_af_dict[af_key]
                                   tumour_af = [i for i in tumour_afs if i[0] == alt][0].split(":")
                               except KeyError:
                                   tumour_af = ["N", "0", "0", "0", "0"]
                                   #print "tumour KEYERROR! %s\n" % tumour_af
                               t_totC, t_refC, t_altC, t_af = tumour_af[1], tumour_af[2], tumour_af[3], tumour_af[4]
                               t_alt = tumour_af[0]
                               writer.writerow([chr, start, ref, alt, n_totC, n_refC, n_altC, n_af, t_alt, t_totC, t_refC, t_altC, t_af])
            else:
                n_af_dict = dict()
                n_af_dict["no_matched_normal"] = "True"   
                print "%s %s normal bam does not exist, merging not required!" % (patient, type)

    return com_af_files


def detect_cluster_job_status(completeion_file_list):
    completed = False
    for file in completeion_file_list:
        if (os.path.exists(file)):
            completed = True
        else:
            completed = False
            break
    return completed

def detect_cluster_jobs(complete_stamps):
    """ detect if job on cluster finised """
    job_status = False
    print "Waiting for cluster jobs to finish!\n"
    while (not job_status):
        time.sleep(10)
        job_status = detect_cluster_job_status(complete_stamps)
    print "All cluster jobs finished? %s\n" % job_status


def make_pileup_scripts(patient_files):
    """ genearate pileup scripts, which also parses pileup results! """
    pileup_scripts = []
    pileup_complete_stamps = []
    for patient in patient_files:
        for status in patient_files[patient]:
            bam_types = dict()
            snv_pos_file = ".".join([patient, "vcf.snp.pos"])
            DNA_bam = patient_files[patient][status]['DNA_bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            bam_types["DNA"] = DNA_bam
            bam_types["RNA"] = RNA_bam
            for type in bam_types:
                if (bam_types[type] != "NA"): 
                    script = make_pileup_script(patient, status, type, snv_pos_file, bam_types[type])
                    pileup_scripts.append( script )
                    pileup_complete_stamp = "_".join([patient, status, type, "pileup.complete"])
                    pileup_complete_stamps.append( pileup_complete_stamp )
    #print pileup_scripts
    return [pileup_scripts, pileup_complete_stamps]


def make_af_dict(patient_files):
    """ 
    Reading in af files for all patients and making af_dicts:
    {'DNA': [{'Axxxxx_unknown_10_100148176_A_G': ['2,0,2,1.00,2,0,2,1.00,2.0,0.0,2.0,1.00,1.000'],
              'Axxxxx_unknown_X_99657566_G_T': ['4,2,2,0.50,4,2,2,0.50,4.0,2.0,2.0,0.50,1.000'],
              'Axxxxx_unknown_X_99941705_T_C': ['5,0,5,1.00,5,0,5,1.00,5.0,0.0,5.0,1.00,1.000']},
             {}],
     'RNA': [{},
             {'Axxxxx_unknown_10_100148176_A_G': ['2,0,2,1.00'],
              'Axxxxx_unknown_X_99941705_T_C': ['5,0,5,1.00']}]}
    """

    AFcounts_afs_wn = dict()
    af_dicts = dict()
    types = ["DNA", "RNA"]
    for type in types:
        bam_key = "_".join([type, "bam"])
        affs_wn = []
        affs_non = []
        af_dicts[type] = []
        AFcounts_afs_wn[type] = []
        for patient in patient_files:
            # see if DNA or RNA normal BAMs available
            if ("normal" in patient_files[patient]):
                nor_bam = patient_files[patient]["normal"][bam_key]
            else:
                nor_bam = "NA"
            if (nor_bam != "NA" ):
                print "%s has %s normal bam!" % (patient, type)
                for status in patient_files[patient]:
                    #bam = patient_files[patient][status][bam_key]
                    if (not status == "normal"): #and (bam != "NA"):
                        aff  = ".".join([patient, status, type, "af.combined.adjusted.fisher"])
                        affs_wn.append(aff)
                        AFcounts_af  = ".".join([patient, status, type, "pileup.AFcounts.af"])
                        AFcounts_afs_wn[type].append( AFcounts_af )
    
            elif (nor_bam == "NA"):
                print "%s does not have %s normal bam!" % (patient, type)
                for status in patient_files[patient]:
                    bam = patient_files[patient][status][bam_key]
                    #if (bam != "NA"):
                    if (not status == "normal") and (bam != "NA"):
                        aff  = ".".join([patient, status, type, "pileup.AFcounts.af"])
                        affs_non.append( aff )
 
        print "%s no normal af files are: %s" % (type, affs_non)
        print "%s with normal af files are: %s" % (type, affs_wn)

        af_non_dict = make_patient_non_af_dict( affs_non )
        #'PASWLN_Primary_3_50879176_A_T': ['18,14,4,0.22']
    
        af_wn_dict = make_patient_wn_af_dict( affs_wn )
        #'PASWAJ_Relapse_8_116635942_C_T': ['96,55,41,0.43,116,61,55,0.47,116.0,61.0,55.0,0.47,0.579,Unknown']

        af_dicts[type] = [af_wn_dict, af_non_dict]

    return [af_dicts, AFcounts_afs_wn]



def combine_sum_af_anno(pairing_status, snv_sum, snv_sum_tmp, DNA_af_wn_dict, RNA_af_wn_dict,
                            DNA_af_non_dict, RNA_af_non_dict, patient_files, sum_header):
    with open(snv_sum,  'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter='\t')
        writer.writerow(sum_header)
        with open(snv_sum_tmp,  'r') as snv_fh:
            records = csv.DictReader(snv_fh,  delimiter='\t')
            header = records.fieldnames
            if (pairing_status == "unpaired"):
                for line in records:
                    line_content = [line[i] for i in header ]
                    patient_id = line['patient_ID']
                    gene = line['gene']
                    chr = line['chromosome']
                    pos = line['position']
                    ref = line['ref_base']
                    alt = line['alt_base']
                    sl = patient_id.split("_")
                    patient = sl[0]
                    status = '_'.join(sl[1:])
                    DNA_bam = patient_files[patient][status]['DNA_bam']
                    RNA_bam = patient_files[patient][status]['RNA_bam']
                    DNA_tc = patient_files[patient][status]['DNA_tc']
                    RNA_tc = patient_files[patient][status]['RNA_tc']
                    variant = "_".join([patient_id, chr, pos, ref, alt])

                    # DNA af info 
                    try:
                        DNA_af = DNA_af_non_dict[variant][0].split(',')
                    except KeyError:
                        if (DNA_bam == "NA"):
                            DNA_af = ["na", "na", "na", "na"]
                        else:
                            DNA_af = ["0", "0", "0", "0"]

                    # RNA af info 
                    try:
                        RNA_af = RNA_af_non_dict[variant][0].split(',')

                    except KeyError:
                        if (RNA_bam == "NA"):
                            RNA_af = ["na", "na", "na", "na"]
                        else:
                            RNA_af = ["0", "0", "0", "0"]

                    final_content = line_content + DNA_af + [DNA_tc] + RNA_af + [RNA_tc]
                    writer.writerow(final_content)

            elif (pairing_status == "paired"):
                 for line in records:
                     line_content = [line[i] for i in header ]
                     patient_id = line['patient_ID']
                     gene = line['gene']
                     chr = line['chromosome']
                     pos = line['position']
                     ref = line['ref_base']
                     alt = line['alt_base']
                     sl = patient_id.split("_")
                     patient = sl[0]
                     status = '_'.join(sl[1:])
                     n_DNA_bam = patient_files[patient]['normal']['DNA_bam']
                     n_RNA_bam = patient_files[patient]['normal']['RNA_bam']
                     DNA_bam = patient_files[patient][status]['DNA_bam']
                     RNA_bam = patient_files[patient][status]['RNA_bam']
                     DNA_tc = patient_files[patient][status]['DNA_tc']
                     RNA_tc = patient_files[patient][status]['RNA_tc']
                     variant = "_".join([patient_id, chr, pos, ref, alt])
 
                     # DNA af info
                     if (n_DNA_bam != "NA"): 
                         try:
                             DNA_af = DNA_af_wn_dict[variant][0].split(',')
                         except KeyError:
                             if (DNA_bam != "NA"): 
                                 DNA_af = ["0", "0", "0", "0" ,"0", "0", "0", "0", "0", "0", "0", "0", "0"]
                             else:
                                 DNA_af = ["na", "na", "na", "na" ,"na", "na", "na", "na", "na", "na", "na", "na", "na"] 
                     elif (n_DNA_bam == "NA"):
                         try:
                             DNA_af_tmp = DNA_af_non_dict[variant][0].split(',')
                             na = ["na", "na", "na", "na" ]
                             DNA_af = na + DNA_af_tmp + na + ["na"]
                              
                         except KeyError:
                             if (DNA_bam != "NA"): 
                                 DNA_af = ["na", "na", "na", "na" ,"0", "0", "0", "0", "na", "na", "na", "na", "na"]
                             else:
                                 DNA_af = ["na", "na", "na", "na" ,"na", "na", "na", "na", "na", "na", "na", "na", "na"]
                     # RNA af info 
                     if (n_RNA_bam != "NA"):
                         try:
                             RNA_af = RNA_af_wn_dict[variant][0].split(',')
                         except KeyError:
                             RNA_af = ["0", "0", "0", "0" ,"0", "0", "0", "0", "0", "0", "0", "0", "0"]
                     elif (n_RNA_bam == "NA"):
                         try:
                             RNA_af_tmp = RNA_af_non_dict[variant][0].split(',')
                             na = ["na", "na", "na", "na" ]
                             RNA_af = na + RNA_af_tmp + na + ["na"]

                         except KeyError:
                             if (RNA_bam != "NA"):
                                 RNA_af = ["na", "na", "na", "na" ,"0", "0", "0", "0", "na", "na", "na", "na", "na"]
                             else:
                                 RNA_af = ["na", "na", "na", "na" ,"na", "na", "na", "na", "na", "na", "na", "na", "na"]
                     final_content = line_content + DNA_af + [DNA_tc] + RNA_af + [RNA_tc]
                     writer.writerow(final_content)


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

def check_file_permission(patient_files):
    strelka_snv_vcfs = []
    strelka_indel_vcfs = []
    DNA_bams = []
    paired_mpileup_vcfs = []
    mutseq_snv_vcfs = []
    DNA_single_vcfs = []
    RNA_single_vcfs = []
    RNA_bams = []
    for patient in patient_files:
        for status in patient_files[patient]:
            strelka_snv_vcf = patient_files[patient][status]['strelka_snv_vcf']
            strelka_indel_vcf = patient_files[patient][status]['strelka_indel_vcf']
            DNA_single_vcf = patient_files[patient][status]['DNA_single_vcf']
            paired_mpileup_vcf = patient_files[patient][status]['paired_mpileup_vcf']
            mutseq_snv_vcf = patient_files[patient][status]['mutseq_snv_vcf']
            DNA_bam = patient_files[patient][status]['DNA_bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            RNA_single_vcf = patient_files[patient][status]['RNA_single_vcf']

            if (strelka_snv_vcf != "NA"):
                strelka_snv_vcfs.append(strelka_snv_vcf)
            if (strelka_indel_vcf != "NA"):
                strelka_indel_vcfs.append(strelka_indel_vcf)
            if (DNA_single_vcf != "NA"):
                DNA_single_vcfs.append(DNA_single_vcf)
            if (RNA_single_vcf != "NA"):
                RNA_single_vcfs.append(RNA_single_vcf)
            if (paired_mpileup_vcf != "NA"):
                paired_mpileup_vcfs.append(paired_mpileup_vcf)
            if (mutseq_snv_vcf != "NA"):
                mutseq_snv_vcfs.append(mutseq_snv_vcf)
            if (DNA_bam != "NA"):
                DNA_bams.append(DNA_bam)
            if (RNA_bam != "NA"):
                RNA_bams.append(RNA_bam)

    print "Checking strelka snv vcf fileis permissions!"
    check_files(strelka_snv_vcfs)
    print "Checking strelka indel vcf fileis permissions!"
    check_files(strelka_indel_vcfs)

    print "Checking paired_mpileup vcf file permissions!"
    check_files(paired_mpileup_vcfs)

    print "Checking single_mpileup vcf file permissions!"
    check_files(DNA_single_vcfs)
    check_files(RNA_single_vcfs)

    print "Checking bam file permissions!"
    check_files(DNA_bams)

    print "Checking mutseq snv vcf file permissions!"
    check_files(mutseq_snv_vcfs)

    print "Checking RNA_bam file permissions!"
    check_files(RNA_bams)



def filter_variants(variant_summary):
    filtered_summary = ".".join([variant_summary, "filtered" ])
    somatic_summary = ".".join([filtered_summary, "somatic" ])
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
                     DNA_n_cov = int(line["n_DNA_cov"])
                     DNA_n_refC = int(line["n_DNA_RefC"])
                     DNA_n_altC = int(line["n_DNA_AltC"])
                     DNA_n_af = float(line["n_DNA_AF"])
                     DNA_t_cov = int(line["t_DNA_cov"])
                     DNA_t_refC = int(line["t_DNA_RefC"])
                     DNA_t_altC = int(line["t_DNA_AltC"])
                     DNA_t_af = float(line["t_DNA_AF"])
                     RNA_n_cov = int(line["n_RNA_cov"])
                     RNA_n_refC = int(line["n_RNA_RefC"])
                     RNA_n_altC = int(line["n_RNA_AltC"])
                     RNA_n_af = float(line["n_RNA_AF"])
                     RNA_t_cov = int(line["t_RNA_cov"])
                     RNA_t_refC = int(line["t_RNA_RefC"])
                     RNA_t_altC = int(line["t_RNA_AltC"])
                     RNA_t_af = float(line["t_RNA_AF"])
                     RNA_t_altref_total = RNA_t_refC + RNA_t_altC
                     content = [line[i] for i in headers]
                     #quality filtering
                     if ((DNA_t_cov > 10 and DNA_t_altC > 3 and DNA_t_af >= 0.1) or
                        (RNA_t_cov > 10 and RNA_t_altC > 3 and RNA_t_af >= 0.1 and RNA_t_altref_total > 0.8*RNA_t_cov)):
                         writer.writerow(content)
                         # somatic filters
                         if ((DNA_n_af < 0.03) or (DNA_n_altC < 2)) and ((RNA_n_af < 0.03) or (RNA_n_altC < 2)):
                             #print "somatic", DNA_n_af, DNA_n_altC, RNA_n_af, RNA_n_altC
                             writer2.writerow(content)
    return [filtered_summary, somatic_summary]


def make_expands_input_file(filtered_summary, patient_files_wn):
    """ 
    d = gene_variant_patients dictionary 
    tc_d = tumor content dictionary
    """
    d = dict()
    tc_d = dict()
    snv_header = ["chr", "startpos", "T_refC", "T_altC","AF_Tumor", "PN_B", "gene"]
    cn_header = ["chr", "startpos", "endpos", "CN_Estimate"]
    # header = ["mutation_id", "ref_counts", "var_counts","normal_cn", "minor_cn", "major_cn", "variant_freq"]
    with open (filtered_summary, 'r')  as handle:
         records = csv.DictReader(handle,  delimiter='\t')
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
             pat = line["patient_ID"]
             if chr.isdigit() and (DNA_n_af < DNA_t_af): # tumour af has to be > n_af
                 if (DNA_n_af < 0.03 or DNA_n_altC < 2):
                     alt_ploidy = 0 # somatic 
                 else:
                     alt_ploidy = 1 # germline
                 variant = ":".join([chr, pos,  str(DNA_t_refC), str(DNA_t_altC), str(DNA_t_af), str(alt_ploidy), gene])
                 try:
                     d[pat].append( variant )
                 except KeyError:
                     #print "key error!"
                     d[pat] = [variant]
    print d

    print "Making expands tsv input files!\n"
    for pat in d:
        snv_file = "".join([pat, ".expands.snv" ])
        with open (snv_file,  'wb') as fh:
            writer = csv.writer( fh, delimiter='\t' )
            writer.writerow( snv_header )
            for mutation in list(set(d[pat])):
                sl =  mutation.split(":")
                writer.writerow(sl)

    print "Generating expands R scripts!\n" 
    wkdir = os.getcwd() 
    wkdir = wkdir + "/"
    for patient in patient_files_wn:
        for status in patient_files_wn[patient]:
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
 
                #copy cn file and add header, cn needs to be absolute copy number
                cnv = patient_files_wn[patient][status]["cnv"]
                shutil.copyfile(cnv, cn_file)
                for line in fileinput.input(cn_file, inplace=True):
                    if fileinput.isfirstline():
                        print '\t'.join(cn_header)
                    print line,
     
def __main__():
    print "Gene summary scripts starts at: %s\n" % datetime.datetime.now()     
    print "-------------------------------------------------------------------------------"
    print "This pipeline does the following:\n \
           -> Combine DNA and RNA variants from mpileup, mutseq, strelka pipeline etc,\n \
           -> Summarize variants by gene and variant level,\n \
           -> Rerun samtools mpileup at default setting for all variant positions, \n \
           -> Parse mpileup results to get base count info, \n \
           -> Combining base count info for tumour-normal pairs if applicable, \n \
           -> Adjust tumour base count according to tumour content, \n \
           -> Perform Fisher Exact test to determine p value, \n \
           -> Annotate variant in summary with base counts and cosmic, \n \
           -> Indicate the caller(s) for each variant.\n"
    print "-------------------------------------------------------------------------------"
            

    parser = argparse.ArgumentParser(description='Summarize variants at both gene and variant level')
    parser.add_argument('-dt','--data_type', help='specify data type as wgs or spc', required=True)
    parser.add_argument('-i','--input_file', help='specify input file', required=True)
    args = vars(parser.parse_args())

    # variant_input_file = "bam_vcf_cnv_path.txt" 
    variant_input_file = args['input_file']


    # impact_types as annotated by snpeff
    impact_types = ["HIGH", "MODERATE"]

    # use this R path when running qsub or locally
    Rscript_path = "/gsc/software/linux-x86_64-centos6/R-3.1.1/bin"
    r_script = "fisherExactTest.r"

    # snv file names for patients with normal
    hm_snv_sum_wn = "high_moderate_SNV_summary_with_normal.txt"
    hm_snv_sum_wn_tmp = ".".join([hm_snv_sum_wn, "tmp"])

    # snv file names for patients without normal
    hm_snv_sum_non = "high_moderate_SNV_summary_no_normal.txt"
    hm_snv_sum_non_tmp = ".".join([hm_snv_sum_non, "tmp"])

    # indel file names for patients with normal
    hm_indel_sum_wn = "high_moderate_INDEL_summary_with_normal.txt"

    # indel file names for patients without normal
    hm_indel_sum_non = "high_moderate_INDEL_summary_no_normal.txt"

    sum_snv_header = ["gene",  "num_patients_gene_level", "num_SNVs_gene_level",
                      "chromosome", "position", "ref_base", "alt_base", 
                      "num_patients_SNV_level", "patient_ID", "snp_ID", "gmaf", 
                      "cosmic64_ID", "snpeff_details", "in_DNA_single_mpileup", "in_paired_mpileup",
                      "in_mutseq", "in_strelka", "in_RNA_single_mpileup"]

    sum_snv_header_non = sum_snv_header + ["t_DNA_cov", "t_DNA_RefC", "t_DNA_AltC", "t_DNA_AF", "DNA_tc",
                                           "t_RNA_cov", "t_RNA_RefC", "t_RNA_AltC", "t_RNA_AF", "RNA_tc"]

    sum_snv_header_wn = sum_snv_header + ["n_DNA_cov", "n_DNA_RefC", "n_DNA_AltC", "n_DNA_AF", 
                                          "t_DNA_cov", "t_DNA_RefC", "t_DNA_AltC", "t_DNA_AF", 
                                          "adj_t_DNA_cov", "adj_t_DNA_RefC", "adj_t_DNA_AltC", "adj_t_DNA_AF",
                                          "DNA_fisher_pvalue", "DNA_tc", 
                                          "n_RNA_cov", "n_RNA_RefC", "n_RNA_AltC", "n_RNA_AF", 
                                          "t_RNA_cov", "t_RNA_RefC", "t_RNA_AltC", "t_RNA_AF", 
                                          "adj_t_RNA_cov", "adj_t_RNA_RefC", "adj_t_RNA_AltC", "adj_t_RNA_AF",
                                          "RNA_fisher_pvalue", "RNA_tc"]

    sum_indel_header = ["gene",  "num_patients_gene_level", "num_INDELs_gene_level", 
                        "chromosome", "position", "ref_base", "alt_base", 
                        "num_patients_INDEL_level", "patient_ID", "snp_ID", "gmaf", 
                        "cosmic64_ID", "snpeff_details", "in_mpileup", "in_strelka", "pileup_cov",
                        "pileup_RefC", "pileup_AltC", "pileup_AF",
                        "strelka_n_Cov", "strelka_n_RefC", "strelka_n_AltC", "strelka_n_AF",
                        "strelka_t_Cov", "strelka_t_RefC", "strelka_t_AltC", "strelka_t_AF"]


    print "Removing completion stamps so that the pipeline can be initiated!\n"
    extension = ['intersect.complete', 'pileup.complete']
    delete_files_by_extension(extension)

    print "Variant input file is:\n%s\n" % (variant_input_file)

    print "-------------------------------------------------------------------------------"
    print "Generating patient_files dictionary!\n" 
    print "-------------------------------------------------------------------------------"
    patient_files = make_patient_files_dict(variant_input_file)
    out_dict = group_patients(patient_files)
    patient_files_wn = out_dict[0]
    patient_files_non = out_dict[1]
    


    print "Patients with matched normal are:\n"
    pprint(patient_files_wn)
    
    print "Patients without matched normal are:\n"
    pprint(patient_files_non)

    print "-------------------------------------------------------------------------------"
    print "Checking bam, vcf file permissions!"
    print "-------------------------------------------------------------------------------"
    check_file_permission(patient_files)


    if args['data_type'] == 'wgs':
         
        print "-------------------------------------------------------------------------------"
        print "Summarizing indels in vcfs for tumours with matching normals!\n"
        print "-------------------------------------------------------------------------------"
        summarize_indels(patient_files_wn, sum_indel_header, hm_indel_sum_wn, impact_types)

        print "-------------------------------------------------------------------------------"
        print "Summarizing indels in vcfs for tumours without normals!\n"
        print "-------------------------------------------------------------------------------"
        summarize_indels(patient_files_non, sum_indel_header, hm_indel_sum_non, impact_types)
        
       
        print "-------------------------------------------------------------------------------"
        print "Summarizing snvs in vcfs for tumours with matching normals!"
        print "-------------------------------------------------------------------------------"
        sum_out_wn = summarize_snvs(patient_files_wn, sum_snv_header, hm_snv_sum_wn_tmp, impact_types)
        patient_snv_wn = sum_out_wn[0]
        pos_files_wn = sum_out_wn[1]
         
        print "-------------------------------------------------------------------------------"
        print "Summarizing snvs in vcfs for tumours without matching_normal!"
        print "-------------------------------------------------------------------------------"
        sum_out_non = summarize_snvs(patient_files_non, sum_snv_header, hm_snv_sum_non_tmp, impact_types)
        patient_snv_non = sum_out_non[0]
        pos_files_non = sum_out_non[1]
         
        print "-------------------------------------------------------------------------------"
        print "Generating mpileup scripts!\n"
        print "-------------------------------------------------------------------------------"
        out = make_pileup_scripts(patient_files) 
        pileup_scripts = out[0]
        pileup_complete_stamps = out[1]
        print pileup_scripts
        #print pileup_complete_stamps

        print "-------------------------------------------------------------------------------"
        print "Qsub mpileup scripts!\n"
        print "-------------------------------------------------------------------------------"
        #qsub_scripts(pileup_scripts)
    
        
        print "-------------------------------------------------------------------------------"
        print "Detecting if pileup jobs on cluster finised!\n"
        print "-------------------------------------------------------------------------------"
        #detect_cluster_jobs(pileup_complete_stamps)
                
        print "Deleting files! \n"
        #delete_files(pileup_scripts)
        #delete_files (pos_files_non)
        #delete_files (pos_files_wn)
        
        print "-------------------------------------------------------------------------------"
        print "Parsing pileup output to get base counts!\n"
        print "-------------------------------------------------------------------------------"
        com_af_files = parse_pileup_output(patient_files, patient_snv_wn, patient_snv_non)
        print "com_af_files are:"
        pprint(com_af_files)
 
        print "-------------------------------------------------------------------------------"
        print "Adjusting tumour base count by tumour content!\n"
        print "-------------------------------------------------------------------------------"
        for type in com_af_files:
            adjust_allele_frequency(type, com_af_files[type], patient_files)

        print "-------------------------------------------------------------------------------"
        print "Performing Fisher Exact Test for tumour/normal pairs!\n"
        print "-------------------------------------------------------------------------------"
        run_fisher_exact_test(Rscript_path, r_script)
        
        
        print "Deleting intermediate files!" 
        extension=['combined.adjusted', 'combined', '.pileup.log']
        #delete_files_by_extension(extension)
        
        print "-------------------------------------------------------------------------------"
        print "Figuring out if DNA and/or RNA bams available!\n"
        print "-------------------------------------------------------------------------------"

        print "-------------------------------------------------------------------------------"
        print "Reading in DNA and RNA af files and making af dictionary! \n"
        print "-------------------------------------------------------------------------------"
        af_out = make_af_dict(patient_files)
        DNA_af_wn_dict = af_out[0]["DNA"][0]
        DNA_af_non_dict = af_out[0]["DNA"][1]
        RNA_af_wn_dict = af_out[0]["RNA"][0]
        RNA_af_non_dict = af_out[0]["RNA"][1]
        DNA_AFcounts_afs_wn = af_out[1]["DNA"]
        RNA_AFcounts_afs_wn = af_out[1]["RNA"]

        print "Deleting AFcounts.af files!\n"
        #delete_files(AFcounts_afs_wn)
 
        print "-------------------------------------------------------------------------------"
        print "Combining snv_non summary, af, and annotation results!\n"
        print "-------------------------------------------------------------------------------"
        pairing_status = 'unpaired'
        combine_sum_af_anno(pairing_status, hm_snv_sum_non, hm_snv_sum_non_tmp, DNA_af_wn_dict, RNA_af_wn_dict,
                            DNA_af_non_dict, RNA_af_non_dict, patient_files, sum_snv_header_non)
    
        
        print "-------------------------------------------------------------------------------"
        print "Combining snv_wn summary, af, and annotation results!\n"
        print "-------------------------------------------------------------------------------"
        pairing_status = 'paired'
        combine_sum_af_anno(pairing_status, hm_snv_sum_wn, hm_snv_sum_wn_tmp, DNA_af_wn_dict, RNA_af_wn_dict,
                            DNA_af_non_dict, RNA_af_non_dict, patient_files, sum_snv_header_wn)
        
        print "Deleting intermediate files!"
        extension=['adjusted.fisher', '.vcf', 'AFcounts.af']
        #delete_files_by_extension(extension)
        
    print "Summarization scripts finished on: %s\n" % datetime.datetime.now()    

if __name__ == '__main__':
    __main__()

