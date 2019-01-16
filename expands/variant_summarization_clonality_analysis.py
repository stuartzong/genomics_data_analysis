import os
import os.path
import time
import datetime
import subprocess
import re
import sys
import glob
import argparse
import csv
from collections import defaultdict
from pprint import pprint
from itertools import islice
import fileinput
import shutil

"""
testing script: expands pipeline (v1.0) for production
To run this script, type the following command:
python  variant_summarization_clonality_analysis.py -dt wgs or spc > summary_log_file.txt

"""
def make_patient_files_dict_new(bam_vcf_files):
    """
    dictionary holds all files: patient -> status -> file_identifier -> file_path
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
            paired_mpileup_vcf = line['paired_mpileup_vcf']
            mutseq_snv_vcf = line['mutseq_snv_vcf']
            bam = line['bam_path']
            strelka_snv_vcf = line['strelka_snv_vcf']
            strelka_indel_vcf = line['strelka_indel_vcf']
            other_vcf = line['other_vcf']
            cnv = line['cnv']
            RNA_bam = line['RNA_bam']
            RNA_single_vcf = line['RNA_single_vcf']
            tmp_dict = {'lib':lib,
                    'single_vcf':single_vcf,
                    'paired_mpileup_vcf': paired_mpileup_vcf,
                    'mutseq_snv_vcf': mutseq_snv_vcf,
                    'bam':bam,
                    'strelka_snv_vcf':strelka_snv_vcf,
                    'strelka_indel_vcf':strelka_indel_vcf,
                    'cnv':cnv,
                    'tc':tc,
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
        #pprint(patient_files)
        return patient_files


def take(n, iterable):
    " Return first n items of the iterable as a list! "
    return list(islice(iterable, n))

def is_other_readable(file_path):
    m = os.stat(file_path).st_mode
    #otherExec  = bool(m & 0001)
    #otherWrite = bool(m & 0002)
    read_status  = bool(m & 0004)
    return read_status
 
def exist(file_path):
    exist_status = os.path.exists(file_path)
    return exist_status
 
def check_files(files):
    missing_files = []
    for file in files:
        print "Checking %s\n" % file
        exist_status = False
        read_status = False
        exist_status = exist(file)
        if (exist_status):
            read_status = is_other_readable(file)
        if (not(exist_status and read_status)):
            missing_files.append(file)
        else:
            print "File exists and readable!\n"
    if not missing_files:
        print "All files exist and are readable!\n" 
    else:
        print "ERROR: The following files are either missing or no read permission!\n"
        print "ATTENTION: Please check to make sure STRELKA has completed for this sample!\n"
        for file in missing_files:
            print file

def make_intersect_script(vcf_file, patient, status, target_region, pipeline):
    """
    Add three extra column to vcf files and 
    Removed header for bedtools interscetBed
    Generate intersect scripts to remove off-target calls
    """
    bed_path = "/home/rcorbett/aligners/bedtools/BEDTools-Version-2.15.0/bin/"
    sam_path = "/home/rcorbett/aligners/samtools/samtools-0.1.17/"
    no_header = ".".join([patient, status, pipeline, "noheader"])
    wkdir = os.getcwd()
    # sometimes no strelka file because of no somatic calls
    exist_status = exist(vcf_file)
    if (exist_status):
        file_size = str(os.stat(vcf_file)).split(',')[6].split('=')[1]
        print no_header, file_size, type(file_size)
        if (int(file_size) != 0):
            try:
                with open(vcf_file, 'r') as handle:
                    with open(no_header, 'w') as writer:
                        for line in handle:
                            if not line.startswith("#"):
                                split = line.strip().split("\t")
                                writer.write("\t".join([split[0], split[1], split[1], line]))
            except IOError:
                print "ERROR: %s exists but can not be read correctly. Network issue?\n" % vcf_file
        else:
            print "ATTENTION: Assume no Strelka somatic calls and continue summarization pipeline.\n"
            with open(no_header, 'wb') as no_writer:
                no_writer.write("100\t1\t1\n")
    
        
    else:
        print "ERROR: %s does not exist. likely due to strelka filters, no calls = no file.\n" % vcf_file
        print "ATTENTION: Assume no Strelka somatic calls and continue summarization pipeline.\n"
        with open(no_header, 'wb') as no_writer:
            no_writer.write("100\t1\t1\n")

    # For some reason, subprocess ignores some lines at end of the file, depreciated this approach.
    # instead generate bash script to intersect

    intersected_file = ".".join([patient, status, pipeline, "intersected"])
    final_isected  = ".".join([patient, status, pipeline, "vcf"])
    bash_file = ".".join([patient, status, pipeline, "intersect.sh"])
    with open(bash_file,  'wb') as writer:
        writer.write("#! /bin/bash\n")
        writer.write("#$ -S /bin/bash\n")
        writer.write("#$ -N bed_intersect\n")
        writer.write("#$ -q transabyss.q\n")
        writer.write("#$ -l mem_token=2G,mem_free=2G,h_vmem=2G\n")
        writer.write("#$ -V\n")
        writer.write("%s%s%s\n" % ("export PATH=", sam_path, ":$PATH"))
        writer.write("%s%s%s\n" % ("export PATH=", bed_path, ":$PATH"))
        writer.write("%s%s\n" % ("cd ", wkdir))

        writer.write("%s %s %s %s %s %s\n" % ("intersectBed -a", no_header, 
                                               "-b", target_region, "-wa >", 
                                               intersected_file))
        writer.write("if [ $? -eq 0 ]\n")
        writer.write("    then\n")
        writer.write("%s %s_%s_%s_%s\n" % ("touch ", patient, status, pipeline,
                     "intersect.complete"))
        writer.write("echo \"intersectBED job finished successfully at:\" `date`\n")
        writer.write("else\n")
        writer.write("    echo \"ERROR: bedtools intersectBED did not finish correctly!\" >>\"summary_log_file.txt\"\n")
        writer.write("".join(["    echo \"Failed patient is: ", patient, "_", status, "\" >>\"summary_log_file.txt\"", "\n"])) 
        writer.write("    exit\n")
        writer.write("fi\n")
        

        writer.write("%s %s %s %s\n" % ("cut -f 4-", intersected_file, ">", 
                                        final_isected))
        writer.write("%s %s\n" % ("rm", intersected_file)) 
        writer.write("%s %s\n" % ("rm", no_header)) 
        writer.write("echo \"intersectBED and post-processing job finished successfully at:\" `date`\n")
    return bash_file


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

def hold_snv_info(transcripts, impacts, patient_status, snv, snp_id, cosmic_id, gmaf, caller, gene_pats, gene_variants):
    selected_transcripts = [k for k in transcripts if any(x in k for x in impacts)]
    # pick the first transcript to get the gene name
    trans_details = selected_transcripts[0]
    gene = trans_details.split("|")[5]
    if (gene == ""):
        gene = "INTERGENIC" 
    # make gene -> patient dict to count how many 
    # patients have variant in this gene
    try:
       gene_pats[gene].append(patient_status)
    except:
       gene_pats[gene] = [patient_status]

    details = ":".join([snp_id, cosmic_id, gmaf, trans_details, caller])

    # put all relevant info into a triple nested dictionary: 
    # gene -> snv -> patient -> details
    try:
       gene_variants[gene][snv][patient_status].append(details)
    except:
       gene_variants[gene][snv][patient_status] = [details]
    return [gene_pats, gene_variants]

def hold_indel_info(transcripts, impacts, patient_status, snv, snp_id, cosmic_id, gmaf, caller, gene_pats, gene_variants, cov_info):
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
       gene_variants[gene][snv][patient_status].append(details)
    except:
       gene_variants[gene][snv][patient_status] = [details]
    return [gene_pats, gene_variants]

def parse_vcf(vcf_file, caller, impacts, gene_pats, gene_variants, snv_pos, snv_list, patient_status):
    print "Parsing: %s, Variant caller is: %s!\n: " %  (vcf_file, caller)
    
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
    

def summarize_snvs(patient_files, sum_header_wn_tmp, hm_sum_wn_tmp, low_sum_wn_tmp):
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
    low_variants_dict = tree()
    hm_gene_dict = dict()
    low_gene_dict = dict()
    patient_snv = dict()
    patient_callers = dict()
    snv_pos_files = []
    hm_impacts = ["HIGH", "MODERATE"]
    low_impacts = ["LOW"]
    for patient in patient_files:
        print "Start parsing all vcf files related to patient: %s\n" % patient
        snv_pos_file = ".".join([patient, "vcf.snp.pos"]) 
        snv_pos_files.append(snv_pos_file)
        snv_pos = []        
        snv_list = []
        for status in patient_files[patient]:
            if ("normal" not in status):
                print "Parsing all vcf files related to patient and status: %s, %s\n" % (patient, status)
                strelka_snv_vcf = patient_files[patient][status]['strelka_snv_vcf']
                single_vcf = patient_files[patient][status]['single_vcf']
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
    
                if (single_vcf != "NA"):
                    caller = "single_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_vcf(single_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, snv_pos, snv_list, patient_status)
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
    #write_summary(low_variants_dict, low_gene_dict, low_sum_wn_tmp, sum_header_wn_tmp, caller1, caller2)
    #write_summary(modifier_variants_dict, modifier_gene_dict, modifier_sum_wn_tmp, sum_header_wn_tmp, caller1, caller2)
    #return [hm_variants_dict, patient_snv, snv_pos_files, low_variants_dict, modifier_variants_dict]
    return [patient_snv, snv_pos_files]
    
def write_summary(variants_dict, gene_dict, summary_file, sum_header_wn_tmp, patient_callers, patient_files):
    #write info into summary file
    writer = open(summary_file, 'w')
    writer.write("\t".join(sum_header_wn_tmp))
    writer.write("\n")
    print "Start writing variants_dict into a summary file: %s\n" % summary_file
    for gene in variants_dict:
        short_patient = [i.split("_")[0] for i in gene_dict[gene]]
        numPatientGeneLevel = len(list(set(short_patient)))
        numSNVs = len(variants_dict[gene].keys())
        for snv in variants_dict[gene]:
            short_numPatientSnvLevel = [i.split("_")[0] for i in variants_dict[gene][snv]]
            numPatientSnvLevel = len(list(set(short_numPatientSnvLevel)))
            snv_details = []
            combination_in_snv = [i for i in variants_dict[gene][snv]]
            #print "short_patients in patient_files: %s\n" % patient_file
            #print "short_numPatientSnvLevel are: %s\n" % short_numPatientSnvLevel
            # all patient status combinations
            patients_snv_level = list(set(short_numPatientSnvLevel))
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
                #print "patient and status_numPatientSnvLevel is %s, %s\n" % (patient, status_numPatientSnvLevel)
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
                in_single_mpileup = "not_in_single_mpileup"
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

                if ("single_mpileup" in patient_callers[patient]):
                    if ("single_mpileup" in  callers):
                        in_single_mpileup = "in_single_mpileup"
                else:
                    in_single_mpileup = "single_mpileup_not_run"

                if ("RNA_single_mpileup" in patient_callers[patient]):
                    if ("RNA_single_mpileup" in  callers):
                        in_RNA_single_mpileup = "in_RNA_single_mpileup"
                else:
                    in_single_mpileup = "RNA_single_mpileup_not_run"

                if ("mutseq" in patient_callers[patient]):
                    if ("mutseq" in callers):
                        in_mutseq = "in_mutseq"
                else:
                    in_mutseq = "mutseq_not_run"

                content = "\t".join([gene, str(numPatientGeneLevel), str(numSNVs),
                                        chr, start, ref, alt, str(numPatientSnvLevel), 
                                        patient, rs_id, gmaf, cosmic_id, snp_eff, in_single_mpileup,
                                        in_paired_mpileup, in_mutseq, in_strelka, in_RNA_single_mpileup])
                writer.write(content)
                writer.write("\n")
            for combination in combination_not_called:    
                patient = combination
                # see if variant called by both mpileup and strelka
                in_paired_mpileup = "not_in_paired_mpileup"
                in_mutseq = "not_in_mutseq"
                in_strelka = "not_in_strelka"
                in_single_mpileup = "not_in_single_mpileup"
                in_RNA_single_mpileup = "not_in_RNA_single_mpileup"
                if ("strelka" not in patient_callers[patient]):
                    in_strelka = "strelka_not_run"
 
                if ("paired_mpileup" not in patient_callers[patient]):
                    in_paired_mpileup = "paired_mpileup_not_run"
 
                if ("single_mpileup" not in patient_callers[patient]):
                    in_single_mpileup = "single_mpileup_not_run"
 
                if ("RNA_single_mpileup" not in patient_callers[patient]):
                    in_RNA_single_mpileup = "RNA_single_mpileup_not_run"

                if ("mutseq" not in patient_callers[patient]):
                    in_mutseq = "mutseq_not_run"

                content = "\t".join([gene, str(numPatientGeneLevel), str(numSNVs),
                                        chr, start, ref, alt, str(numPatientSnvLevel), 
                                        combination, rs_id, gmaf, cosmic_id, snp_eff, in_single_mpileup,
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


def summarize_indels(patient_files, sum_header_wn_tmp, hm_sum_wn_tmp, low_sum_wn_tmp):
    tree = lambda: defaultdict(tree)
    hm_variants_dict = tree()                                                                                         
    low_variants_dict = tree()
    modifier_variants_dict = tree()
    hm_gene_dict = dict()
    low_gene_dict = dict()
    modifier_gene_dict = dict()
    hm_indel_dict = dict()
    low_indel_dict = dict()
    patient_callers = dict()    
    modifier_indel_dict = dict()
    hm_impacts = ["HIGH", "MODERATE"]
    low_impacts = ["LOW"]

    for patient in patient_files:
        print "Start parsing all vcf files related to patient: %s\n" % patient
        for status in patient_files[patient]:
            if ("normal" not in status):
                print "Parsing all vcf files related to patient and status: %s, %s\n" % (patient, status)
                strelka_indel_vcf = patient_files[patient][status]['strelka_indel_vcf']
                single_vcf = patient_files[patient][status]['single_vcf']
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

                if (single_vcf != "NA"):
                    caller = "single_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_mpileup_indel(single_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, patient_status)
                    hm_gene_dict = two_dicts[0]
                    hm_variants_dict = two_dicts[1]

                if (RNA_single_vcf != "NA"):
                    caller = "RNA_single_mpileup"
                    try:
                        patient_callers[patient_status].append(caller)
                    except KeyError:
                        patient_callers[patient_status] = [caller]
                    two_dicts = parse_mpileup_indel(single_vcf, caller, hm_impacts, hm_gene_dict, hm_variants_dict, patient_status)
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
    print "Parsing: %s, Variant caller is: %s!\n: " %  (vcf_file, caller)
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
    print "Parsing: %s, Variant caller is: %s!\n: " %  (vcf_file, caller)
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
        numPatientGeneLevel = len(list(set(short_patient)))
        numSNVs = len(variants_dict[gene].keys())
        for snv in variants_dict[gene]:
            short_numPatientSnvLevel = [i.split("_")[0] for i in variants_dict[gene][snv]]
            numPatientSnvLevel = len(list(set(short_numPatientSnvLevel)))
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
                in_single_mpileup = "not_in_single_mpileup"
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

                if ("single_mpileup" in patient_callers[patient]):
                    if ("single_mpileup" in  callers):
                        in_single_mpileup = "in_single_mpileup"
                else:
                    in_single_mpileup = "single_mpileup_not_run"

                if ("RNA_single_mpileup" in patient_callers[patient]):
                    if ("RNA_single_mpileup" in  callers):
                        in_RNA_single_mpileup = "in_RNA_single_mpileup"
                else:
                    in_RNA_single_mpileup = "RNA_single_mpileup_not_run"

                content = "\t".join([gene, str(numPatientGeneLevel), str(numSNVs),
                                        chr, start, ref, alt, str(numPatientSnvLevel), 
                                        patient, rs_id, gmaf, cosmic_id, snp_eff, in_single_mpileup,
                                        in_paired_mpileup, in_strelka, in_RNA_single_mpileup, cov_info])
                writer.write(content)
                writer.write("\n")
  

def adjust_allele_frequency(com_af_files,patient_files):
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
                tc_tmp = patient_files[patient][status]['tc']
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

  
def make_gene_full_name_dict(tabFile):
    """tabFie = /home/szong/projects/resource/ens69/gene_info_ens69_nina.txt. To add full gene name"""
    outDict=dict()
    with open (tabFile) as handle:
        for line in handle:
           sl = line.rstrip().split('\t')
           geneSymbol = sl[2]
           geneFullName = [sl[4]]
           #if geneSymbol in outDict:
           try:
             outDict[geneSymbol].append(geneFullName)
           #if geneSymbol not in outDict:
           #else:
           except KeyError:
             outDict[geneSymbol] = geneFullName
        #print outDict
        handle.close()
        return outDict

def make_cosmic_annotation_dict(infile):
   """infile = /home/szong/projects/resource/cosmic/CosmicMutantExport.v71.tsv,  To add COSMIC annotations"""
   # this sometime does not work properly because the way we report base change is different from what is 
   #in the cosmic database v71. for example for cosm521: we report C>T change , it is G>A in cosmic v71
   outDict = dict()
   with open(infile,  'r') as handle:
      records = csv.DictReader(handle,  delimiter='\t')
      for line in records:
         mutation_tmp = line['Mutation_CDS']
         if ">" in mutation_tmp:
            ref, alt = mutation_tmp[-3], mutation_tmp[-1]
            pos = "_".join(re.split('[:|-]', line['Mutation_GRCh37_genome_position'])[0:2])
            gene = line['Gene_name'].split('_')[0]
            key = "_".join([gene, pos, ref, alt])
            value = ";".join([line['Primary_site'], line['Site_subtype'], 
                             line['Primary_histology'], line['Histology_subtype'], 
                             line['Mutation_Description'], line['SNP'], 
                             line['FATHMM_prediction'], line['Mutation_somatic_status'], 
                             line['Sample_name']])
            try:
               outDict[key] = [value]
            except KeyError:
               print "key not in the dictionary\n"
      handle.close()
      #print "cos_var_annos is %s.", outDict
      return outDict


def make_patient_wn_af_dict(af_files):
   patient_af_dict = dict()
   for file in af_files:
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


def remove_offtargets(patient_files, target_region):    
    #remove vcf header and run intersectBed to remove off target regions
    #generate dict to hold normal and no-normal vcf files
    intersect_scripts = []
    complete_stamp = []
    patient_wn_snv_vcfs = dict()
    patient_wn_indel_vcfs = dict()
    patient_non_vcfs = dict()
    for patient in patient_files:
        if ("normal" in patient_files[patient]):
            print "patient has normal: ", patient, patient_files[patient].keys()
            for status in patient_files[patient]:
                if (not status == "normal"):
                    pipelines = ["pileup", "strelka_snv", "strelka_indel"]
                    paired_vcf = patient_files[patient][status][2]
                    strelka_snv_vcf = patient_files[patient][status][3]
                    strelka_indel_vcf = patient_files[patient][status][7]
                    for pipeline in pipelines:
                        stamp_file = "_".join([patient, status, pipeline,
                                               "intersect.complete"])
                        final_inters_file = ".".join([patient, status, pipeline, "vcf"])

                        if ("pileup" in pipeline):
                            vcf = paired_vcf
                            try:
                                patient_wn_snv_vcfs[patient].append(final_inters_file)
                                patient_wn_indel_vcfs[patient].append(final_inters_file)
                            except KeyError:
                                patient_wn_snv_vcfs[patient] = [final_inters_file]
                                patient_wn_indel_vcfs[patient] = [final_inters_file]
                        elif ("strelka_snv" in pipeline):
                            vcf = strelka_snv_vcf
                            try:
                                patient_wn_snv_vcfs[patient].append(final_inters_file)
                            except KeyError:
                                patient_wn_snv_vcfs[patient] = [final_inters_file]
                        elif ("strelka_indel" in pipeline):
                            vcf = strelka_indel_vcf
                            try:
                                patient_wn_indel_vcfs[patient].append(final_inters_file)
                            except KeyError:
                                patient_wn_indel_vcfs[patient] = [final_inters_file]
                        else:
                            print "Error: Invalid vcf file!\n"
                        script = make_intersect_script(vcf, patient, status,
                                                       target_region, pipeline)
                        intersect_scripts.append(script)
                        complete_stamp.append(stamp_file)
    
        else:
            print "patient no normal: ", patient, patient_files[patient].keys()
            pipeline = "pileup"
            for status in patient_files[patient]:
                single_vcf = patient_files[patient][status][1]
                stamp_file = "_".join([patient, status, pipeline,
                                       "intersect.complete"])
                final_inters_file = ".".join([patient, status, pipeline, "vcf"])
                try:
                    patient_non_vcfs[patient].append(final_inters_file)
                except KeyError:
                    patient_non_vcfs[patient] = [final_inters_file]

                script = make_intersect_script(single_vcf, patient, status,
                                               target_region, pipeline)
                intersect_scripts.append(script)
                complete_stamp.append(stamp_file)
    
    ##print intersect_scripts
    #print complete_stamp
    return [intersect_scripts, complete_stamp,
            patient_wn_snv_vcfs, patient_wn_indel_vcfs, patient_non_vcfs]


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
            bam = patient_files[patient][status]['bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            bam_types["DNA"] = bam
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
                    print "For %s %s tumor, example snvs are:" % (patient, status)
                    snp_len = len(snp_list)
                    if (snp_len >= 5):
                        for x in range(2):
                            print snp_list[x]
                        print "\n"
                    else:
                       for x in (snp_list):
                           print x
                       print "\n"
                    # af_dicts hold all af_dict for each patient_status
                    try:
                        #print patient, status
                        af_dicts[type][patient][status] = af_dict
                    except KeyError:
                        #print "key error!"
                        if patient not in af_dicts[type]:
                            af_dicts[type][patient] = {}
                        if status not in af_dicts[type][patient]:
                            af_dicts[type][patient][status] = af_dict
    
        #print "Merging af files for tumour/normal pairs!"
        if ("normal" in patient_files[patient]):
            nor_bams = dict()
            nor_bams["DNA"] = patient_files[patient]["normal"]['bam']
            nor_bams["RNA"] = patient_files[patient]["normal"]['RNA_bam']
            
            for type in nor_bams:    
                if (nor_bams[type] != "NA"):
                    n_af_dict = af_dicts[type][patient]["normal"]
                else:
                    n_af_dict = dict()
                    n_af_dict["no_matched_normal"] = "True"   
                for status in patient_files[patient]:
                    if (not status == "normal"):
                        print "merging %s %s %s with its corresponding normal!" % (patient,status, type)
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
            print "\n"
        else:
            print "%s has no normal and merging is not required!" % patient

    return com_af_files


def merge_normal_tumor_afs(snp_list, n_af_dict, t_af_dict, com_af_file):
    with open(com_af_file,  'wb') as opf:
        writer = csv.writer(opf,  delimiter='\t')
        #print "snp_list is: ", snp_list
        for snp in snp_list:
           sp_snp = snp.split(':')
           chr, start, ref = sp_snp[0], sp_snp[1], sp_snp[2]
           alt = sp_snp[3]
           af_key = ":".join([chr, start, ref])
           #normal DNA counts
           try:
               normal_afs = n_af_dict[af_key]
               normal_af = [i for i in normal_afs if i[0] == alt][0].split(":")
               #print "normal_af is: %s, %s " % (type(normal_af), normal_af)
           except KeyError:
               normal_af = ["N", "0", "0", "0", "0"]
               #print "normal KEYERROR! %s \n" % normal_af
           n_totC, n_refC, n_altC, n_af = normal_af[1], normal_af[2], normal_af[3], normal_af[4]
           # tumor DNA counts
           try:
               tumour_afs = t_af_dict[af_key]
               tumour_af = [i for i in tumour_afs if i[0] == alt][0].split(":")
           except KeyError:
               tumour_af = ["N", "0", "0", "0", "0"]
               #print "tumour KEYERROR! %s\n" % tumour_af
           t_totC, t_refC, t_altC, t_af = tumour_af[1], tumour_af[2], tumour_af[3], tumour_af[4]
           t_alt = tumour_af[0]
           writer.writerow([chr, start, ref, alt, n_totC, n_refC, n_altC, n_af, t_alt, t_totC, t_refC, t_altC, t_af])








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
            bam = patient_files[patient][status]['bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            bam_types["DNA"] = bam
            bam_types["RNA"] = RNA_bam
            for type in bam_types:
                if (bam_types[type] != "NA"): 
                    script = make_pileup_script(patient, status, type, snv_pos_file, bam_types[type])
                    pileup_scripts.append( script )
                    pileup_complete_stamp = "_".join([patient, status, type, "pileup.complete"])
                    pileup_complete_stamps.append( pileup_complete_stamp )
    print pileup_scripts
    return [pileup_scripts, pileup_complete_stamps]


def make_af_dict(patient_files):
    """ Reading in af files and making af dictionary """
    AFcounts_afs_wn = dict()
    af_dicts = dict()
    types = ["DNA", "RNA"]
    for type in types:
        affs_wn = []
        affs_non = []
        af_dicts[type] = []
        AFcounts_afs_wn[type] = []
        for patient in patient_files:
            if ( "normal" in patient_files[patient] ):
                for status in patient_files[patient]:
                    if (not status == "normal"):
                        aff  = ".".join([patient, status, type, "af.combined.adjusted.fisher"])
                        affs_wn.append( aff )
                        AFcounts_af  = ".".join([patient, status, type, "pileup.AFcounts.af"])
                        AFcounts_afs_wn[type].append( AFcounts_af )
    
            elif ( "normal" not in patient_files[patient] ):
                for status in patient_files[patient]:
                    aff  = ".".join([patient, status, type, "pileup.AFcounts.af"])
                    affs_non.append( aff )
        #'PASWLN_Primary_3_50879176_A_T': ['18,14,4,0.22']
        af_non_dict = make_patient_non_af_dict( affs_non )
    
        #'PASWAJ_Relapse_8_116635942_C_T': ['96,55,41,0.43,116,61,55,0.47,116.0,61.0,55.0,0.47,0.579,Unknown']
        af_wn_dict = make_patient_wn_af_dict( affs_wn )
        af_dicts[type] = [af_wn_dict, af_non_dict]
        #print "af_wn_dict is:\n %s\n" % af_wn_dict
    #return [af_wn_dict, af_non_dict, AFcounts_afs_wn]
    return [af_dicts, AFcounts_afs_wn]


def combine_sum_af_anno(pairing_status, snv_sum, snv_sum_tmp, DNA_af_dict, RNA_af_dict, patient_files, sum_header):
    if (pairing_status == "unpaired"):
        with open(snv_sum,  'wb') as csvfile:
            non_writer = csv.writer(csvfile, delimiter='\t')
            non_writer.writerow(sum_header)
            with open(snv_sum_tmp,  'r') as snv_fh:
                 records = csv.DictReader(snv_fh,  delimiter='\t')
                 header = records.fieldnames
                 for line in records:
                     line_content = [line[i] for i in header ]
                     patient = line['patientID']
                     gene = line['gene']
                     chr = line['chromosome']
                     pos = line['position']
                     ref = line['referenceBase']
                     alt = line['alternativeBase']
                     #print "dddddddddd", patient
                     sp_patient = patient.split("_")
                     tar_pat = sp_patient[0]
                     status = "_".join(sp_patient[1:])
                     tc = patient_files[tar_pat][status]['tc']
                     variant = "_".join([patient, chr, pos, ref, alt])
                     # DNA af info 
                     try:
                       DNA_af = DNA_af_dict[variant][0].split(',')
                     except KeyError:
                        DNA_af = ["0", "0", "0", "0"]

                     # RNA af info 
                     try:
                       RNA_af = RNA_af_dict[variant][0].split(',')

                     except KeyError:
                        RNA_af = ["0", "0", "0", "0"]

                     final_content = line_content + DNA_af + RNA_af + [tc]
                     non_writer.writerow(final_content)

    elif (pairing_status == "paired"):
        with open(snv_sum,  'wb') as csvfile:
            writer = csv.writer(csvfile, delimiter='\t')
            writer.writerow(sum_header)
            with open(snv_sum_tmp,  'r') as snv_fh:
                 records = csv.DictReader(snv_fh,  delimiter='\t')
                 header = records.fieldnames
                 for line in records:
                     line_content = [line[i] for i in header ]
                     patient = line['patientID']
                     gene = line['gene']
                     chr = line['chromosome']
                     pos = line['position']
                     ref = line['referenceBase']
                     alt = line['alternativeBase']
                     sp_patient = patient.split("_")
                     #print "dddddddddd", patient
                     #print af_dict
                     tar_pat = sp_patient[0]
                     status = sp_patient[1]
                     tc = patient_files[tar_pat][status]['tc']
                     variant = "_".join([patient, chr, pos, ref, alt])
                     #cosmic_key = "_".join([gene, chr, pos, ref, alt])
    
                     # DNA af info 
                     try:
                       DNA_af = DNA_af_dict[variant][0].split(',')
                     except KeyError:
                        DNA_af = ["0", "0", "0", "0" ,"0", "0", "0", "0", "0", "0", "0", "0", "0"]

                     # RNA af info 
                     try:
                       RNA_af = RNA_af_dict[variant][0].split(',')

                     except KeyError:
                        RNA_af = ["0", "0", "0", "0" ,"0", "0", "0", "0", "0", "0", "0", "0", "0"]
                
                     # cosmic annotation 
                     #try:
                     #   cosmic_anno = cos_var_annos[cosmic_key][0]
                     #except KeyError:
                     #   cosmic_anno = "not_in_cosmic71"
                     #final_content = line_content + af + [tc] + [cosmic_anno]
                     final_content = line_content + DNA_af + RNA_af + [tc]
                     writer.writerow(final_content)



def make_a_list (file):
    any_list = []
    with open(file, 'r') as fh: 
        for line in fh:
           line = line.replace("\n", "")
           any_list.append(line)
    #print any_list
    return any_list
 
 
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
             snp_id = line["SNPID"]
             cosm_id = line["COSMIC64ID"]
             t_refC = line["tumourRefBaseCount"]
             t_altC = line["tumourAltBaseCount"]
             t_tot =line["tumourSequencingCoverage"]
             t_af = line["tumourAlleleFrequency"]
             if (float(t_af) > 0.1 and
                 ("X" not in chr) and ("Y" not in chr) and ("MT" not in chr) and 
                  int(t_tot) >10): # tumour af has to be > n_af
                 if ("not_in_cosmic64" not in cosm_id):
                    alt_ploidy = 0
                 elif ("novel_snp" in snp_id):
                     alt_ploidy = 0 # somatic 
                 else:
                     alt_ploidy = 1 # germline
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
    for patient in patient_files:
        cn_file = "".join([wkdir, "/", patient, ".expands.cn"])
        snv_file = "".join([wkdir, "/", patient, ".expands.snv" ])
        dir ="".join([wkdir, "/", patient])
        r_script = ".".join([ patient, "r" ])
        with open ( r_script,  'wb' ) as info_fh:
            info_fh.write( "".join(["setwd(\"",wkdir,"\")\n"]))
            info_fh.write( "".join(["dir.create(\"", patient, "\", showWarnings = TRUE, recursive = FALSE, mode =""\"0777\")\n"]))
            info_fh.write( "".join(["setwd(\"",dir,"\")\n"]))
            info_fh.write( "library(expands)\n")
            info_fh.write( "\n")
            info_fh.write( "".join(["runExPANdS(\"", snv_file, "\",\"", cn_file, "\",maxScore=2.5, max_PM=6, min_CellFreq=0.1, precision=NA,plotF=2,snvF=\"out.expands\",maxN=8000,region=NA)\n"]) )    
            info_fh.write( "".join(["setwd(\"",wkdir,"\")\n"]))


def copy_cn_file(patient, cn_file):
    """ Copy cn file and add header, cn has to be absolute copy number! """
    cn_header = ["chr", "startpos", "endpos", "CN_Estimate"]
    local_cn_file = ".".join([patient, "expands.cn"])
    shutil.copyfile(cn_file, local_cn_file)
    for line in fileinput.input(local_cn_file, inplace=True):
        if fileinput.isfirstline():
            print '\t'.join(cn_header)
        print line,

def group_patients(patient_files):
    patient_files_wn = dict()
    patient_files_non = dict()
    for patient in patient_files:
        tmp_dict = patient_files[patient]
        if ("normal" in patient_files[patient]):
            patient_files_wn[patient] = tmp_dict
        else:
            patient_files_non[patient] = tmp_dict
    #pprint(patient_files_wn)
    #pprint(patient_files_non)
    return [patient_files_wn, patient_files_non]  

def check_file_permission_availability(patient_files):
    strelka_snv_vcfs = []
    strelka_indel_vcfs = []
    bams = []
    paired_mpileup_vcfs = []
    mutseq_snv_vcfs = []
    single_vcfs = []
    RNA_single_vcfs = []
    RNA_bams = []
    for patient in patient_files:
        for status in patient_files[patient]:
            strelka_snv_vcf = patient_files[patient][status]['strelka_snv_vcf']
            strelka_indel_vcf = patient_files[patient][status]['strelka_indel_vcf']
            single_vcf = patient_files[patient][status]['single_vcf']
            paired_mpileup_vcf = patient_files[patient][status]['paired_mpileup_vcf']
            mutseq_snv_vcf = patient_files[patient][status]['mutseq_snv_vcf']
            bam = patient_files[patient][status]['bam']
            RNA_bam = patient_files[patient][status]['RNA_bam']
            RNA_single_vcf = patient_files[patient][status]['RNA_single_vcf']

            if (strelka_snv_vcf != "NA"):
                strelka_snv_vcfs.append(strelka_snv_vcf)
            if (strelka_indel_vcf != "NA"):
                strelka_indel_vcfs.append(strelka_indel_vcf)
            if (single_vcf != "NA"):
                single_vcfs.append(single_vcf)
            if (RNA_single_vcf != "NA"):
                RNA_single_vcfs.append(RNA_single_vcf)
            if (paired_mpileup_vcf != "NA"):
                paired_mpileup_vcfs.append(paired_mpileup_vcf)
            if (mutseq_snv_vcf != "NA"):
                mutseq_snv_vcfs.append(mutseq_snv_vcf)
            if (bam != "NA"):
                bams.append(bam)
            if (RNA_bam != "NA"):
                RNA_bams.append(RNA_bam)

    print "Checking strelka snv and indel vcf fileis permissions!\n"
    check_files(strelka_snv_vcfs)
    check_files(strelka_indel_vcfs)

    print "Checking paired_mpileup vcf file permissions!\n"
    check_files(paired_mpileup_vcfs)

    print "Checking single_mpileup vcf file permissions!\n"
    check_files(single_vcfs)
    check_files(RNA_single_vcfs)

    print "Checking bam file permissions!\n"
    check_files(bams)

    print "Checking mutseq snv vcf file permissions!\n"
    check_files(mutseq_snv_vcfs)

    print "Checking RNA_bam file permissions!\n"
    check_files(RNA_bams)



def filter_variants(variant_summary):
    filtered_summary = ".".join([variant_summary, "filtered.somatic" ])
    """ d = gene_variant_patients dictionary """
    d = dict()
    with open (filtered_summary,  'wb') as fh:
        writer = csv.writer( fh, delimiter='\t' )
        print "The variant summary file is: %s.\n" % variant_summary
        with open (variant_summary, 'r') as handle:
             records = csv.DictReader(handle,  delimiter='\t')
             headers = records.fieldnames
             writer.writerow(headers)
             for line in records:
                 gene = line['gene']
                 chr = line["chromosome"]
                 pos = line["position"]
                 ref = line["referenceBase"]
                 alt = line["alternativeBase"]
                 DNA_n_cov = int(line["normalDNASequencingCoverage"])
                 DNA_n_refC = int(line["normalDNARefBaseCount"])
                 DNA_n_altC = int(line["normalDNAAltBaseCount"])
                 DNA_n_af = float(line["normalDNAAlleleFrequency"])
                 DNA_t_cov = int(line["tumourDNASequencingCoverage"])
                 DNA_t_refC = int(line["tumourDNARefBaseCount"])
                 DNA_t_altC = int(line["tumourDNAAltBaseCount"])
                 DNA_t_af = float(line["tumourDNAAlleleFrequency"])
                 RNA_n_cov = int(line["normalRNASequencingCoverage"])
                 RNA_n_refC = int(line["normalRNARefBaseCount"])
                 RNA_n_altC = int(line["normalRNAAltBaseCount"])
                 RNA_n_af = float(line["normalRNAAlleleFrequency"])
                 RNA_t_cov = int(line["tumourRNASequencingCoverage"])
                 RNA_t_refC = int(line["tumourRNARefBaseCount"])
                 RNA_t_altC = int(line["tumourRNAAltBaseCount"])
                 RNA_t_af = float(line["tumourRNAAlleleFrequency"])
                 RNA_t_altref_total = RNA_t_refC + RNA_t_altC
                 content = [line[i] for i in headers]
                 #quality filtering
                 if ((DNA_t_cov > 10 and DNA_t_altC > 3 and DNA_t_af >= 0.1) or
                    (RNA_t_cov > 10 and RNA_t_altC > 3 and RNA_t_af >= 0.1 and RNA_t_altref_total > 0.8*RNA_t_cov)):
                     # somatic filters
                     #if (DNA_n_af < 0.03) or (DNA_n_altC < 2) or (RNA_n_af < 0.03) or (RNA_n_altC < 2):
                     if ((DNA_n_af < 0.03) or (DNA_n_altC < 2)) and ((RNA_n_af < 0.03) or (RNA_n_altC < 2)):
                         print "somatic", DNA_n_af, DNA_n_altC, RNA_n_af, RNA_n_altC
                         writer.writerow(content)
    return filtered_summary

def write_list(lst, file):
    with open (file, 'wb') as writer:
        variants = list(set(lst))
        for variant in variants:
            sl = re.split('[:_]', variant)
            writer.write("\t".join(sl))
            writer.write('\n')



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
             ref = line["referenceBase"]
             alt = line["alternativeBase"]
             DNA_n_altC = int(line["normalDNAAltBaseCount"])
             DNA_n_af = float(line["normalDNAAlleleFrequency"])
             DNA_t_refC = int(line["tumourDNARefBaseCount"])
             DNA_t_altC = int(line["tumourDNAAltBaseCount"])
             DNA_t_af = float(line["tumourDNAAlleleFrequency"])
             pat = line["patientID"]
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
    print "Expands scripts starts at: %s\n" % datetime.datetime.now()     
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
            

    parser = argparse.ArgumentParser(description='Summarize variants at both gene and variant level')
    parser.add_argument('-dt','--data_type', help='specify data type as wgs or spc', required=True)
    args = vars(parser.parse_args())

    variant_input_file = "bam_vcf_cnv_path.txt" 
    #input_files = "/projects/trans_scratch/validations/workspace/szong/" \
    #              "AML_capture/all_patients/junk/vcf_bam_all_samples_new_20150924.csv.test"
    cosmic_anno = "/home/szong/projects/resource/cosmic/CosmicMutantExport.v71.tsv"

    # use this R path when running qsub or locally
    Rscript_path = "/gsc/software/linux-x86_64-centos6/R-3.1.1/bin"
    r_script = "fisherExactTest.r"

    # snv file names for patients with normal
    hm_snv_sum_wn = "high_moderate_SNV_summary_with_normal.txt"
    low_snv_sum_wn = "low_SNV_summary_with_normal.txt"
    hm_snv_sum_wn_tmp = ".".join([hm_snv_sum_wn, "tmp"])
    low_snv_sum_wn_tmp = ".".join([low_snv_sum_wn, "tmp"])

    # snv file names for patients without normal
    hm_snv_sum_non = "high_moderate_SNV_summary_no_normal.txt"
    low_snv_sum_non = "low_SNV_summary_no_normal.txt"
    hm_snv_sum_non_tmp = ".".join([hm_snv_sum_non, "tmp"])
    low_snv_sum_non_tmp = ".".join([low_snv_sum_non, "tmp"])

    # indel file names for patients with normal
    hm_indel_sum_wn = "high_moderate_INDEL_summary_with_normal.txt"
    low_indel_sum_wn = "low_INDEL_summary_with_normal.txt"
    modifier_indel_sum_wn = "modifier_INDEL_summary_with_normal.txt"

    # indel file names for patients without normal
    hm_indel_sum_non = "high_moderate_INDEL_summary_no_normal.txt"
    low_indel_sum_non = "low_INDEL_summary_no_normal.txt"
    modifier_indel_sum_non = "modifier_INDEL_summary_no_normal.txt"


    sum_snv_header = ["gene",  "numPatientGeneLevel", "numSNVsForThisGene",
                            "chromosome", "position", "referenceBase", "alternativeBase", 
                            "numPatientSnvLevel", "patientID", "SNPID", "globalMinorAlleleFrequency", 
                            "COSMIC64ID", "impactDetails", "in_single_mpileup", "in_paired_mpileup",
                            "in_mutseq", "in_strelka", "in_RNA_single_mpileup"]
    sum_snv_header_non = sum_snv_header + ["tumourDNASequencingCoverage", "tumourDNARefBaseCount", 
                             "tumourDNAAltBaseCount", "tumourDNAAlleleFrequency",
                             "tumourRNASequencingCoverage", "tumourRNARefBaseCount", 
                             "tumourRNAAltBaseCount", "tumourRNAAlleleFrequency",
                             "tumourContent"]
    sum_snv_header_wn = sum_snv_header + ["normalDNASequencingCoverage", "normalDNARefBaseCount", 
                              "normalDNAAltBaseCount", "normalDNAAlleleFrequency", 
                              "tumourDNASequencingCoverage", "tumourDNARefBaseCount",
                              "tumourDNAAltBaseCount", "tumourDNAAlleleFrequency", 
                              "adjustedDNATumourSequencingCoverage", "adjustedDNATumourRefBaseCount", 
                              "adjustedDNATumourAltBaseCount", "adjustedDNATumourAlleleFrequency",
                              "DNAfisherExactPvalue", 
                              "normalRNASequencingCoverage", "normalRNARefBaseCount", 
                              "normalRNAAltBaseCount", "normalRNAAlleleFrequency", 
                              "tumourRNASequencingCoverage", "tumourRNARefBaseCount",
                              "tumourRNAAltBaseCount", "tumourRNAAlleleFrequency", 
                              "adjustedTumourRNASequencingCoverage", "adjustedTumourRNARefBaseCount", 
                              "adjustedTumourRNAAltBaseCount", "adjustedTumourRNAAlleleFrequency",
                              "RNAfisherExactPvalue", "tumourContent"]
    sum_indel_header = ["gene",  "numPatientGeneLevel", "numINDELsForThisGene", 
                            "chromosome", "position", "referenceBase", "alternativeBase", 
                            "numPatientIndelLevel", "patientID", "SNPID", "globalMinorAlleleFrequency", 
                            "COSMIC64ID", "impactDetails", "in_mpileup", "in_strelka", "pileup_cov",
                            "pileup_refCount", "pileup_altCount", "pileup_alleleFrequency",
                            "strelka_normalCov", "strelka_normalRefCount", "strelka_normalAltCount", "strelka_normalAlleleFrequency",
                            "strelka_tumourCov", "strelka_tumourRefCount", "strelka_tumourAltCount", "strelka_tumourAlleleFrequency"]


    print "Removing completion stamps so that the pipeline can be initiated!\n"
    extension = ['intersect.complete', 'pileup.complete']
    delete_files_by_extension(extension)

    print "Variant input file is:\n%s\n" % (variant_input_file)

    print "Generating patient_files dictionary!\n" 
    patient_files = make_patient_files_dict_new(variant_input_file)
    out_dict = group_patients(patient_files)
    patient_files_wn = out_dict[0]
    patient_files_non = out_dict[1]
    
    print "Patients with matched normal are:\n"
    pprint(patient_files_wn)
    
    print "Patients without matched normal are:\n"
    pprint(patient_files_non)

    print "Making file lists for permission and availability check!\n"
    check_file_permission_availability(patient_files)


    if args['data_type'] == 'wgs':
        '''
        print "Summarizing indels in vcfs for tumours with matching normals!\n"
        summarize_indels(patient_files_wn, sum_indel_header, hm_indel_sum_wn, low_indel_sum_wn)

        print "Summarizing indels in vcfs for tumours without normals!\n"
        summarize_indels(patient_files_non, sum_indel_header, hm_indel_sum_non, low_indel_sum_non)
        
       
        #print "patient_files are: %s\n" % patient_files
        print "Summarizing snvs in vcfs for tumours with matching normals!"
        sum_out_wn = summarize_snvs(patient_files_wn, sum_snv_header, hm_snv_sum_wn_tmp, low_snv_sum_wn_tmp)
        patient_snv_wn = sum_out_wn[0]
        pos_files_wn = sum_out_wn[1]
        print "patient_snv_wn and pos_files_wn are: %s\n%s\n" % (patient_snv_wn ,pos_files_wn)
         
        print "Summarizing snvs in vcfs for tumours without matching_normal!"
        sum_out_non = summarize_snvs(patient_files_non, sum_snv_header, hm_snv_sum_non_tmp, low_snv_sum_non_tmp)
        patient_snv_non = sum_out_non[0]
        pos_files_non = sum_out_non[1]
        #print "patient_snv_non and pos_files_non are: %s\n%s\n" % (patient_snv_non ,pos_files_non)
         
        print "Generating mpileup scripts!\n"
        out = make_pileup_scripts(patient_files) 
        pileup_scripts = out[0]
        pileup_complete_stamps = out[1]
        print pileup_scripts
        print pileup_complete_stamps

        print "Qsub mpileup scripts!\n"
        qsub_scripts(pileup_scripts)
    
        
        print "Detecting if pileup jobs on cluster finised!\n"
        detect_cluster_jobs(pileup_complete_stamps)
                
        print "Deleting files! \n"
        delete_files(pileup_scripts)
        delete_files(pileup_complete_stamps)
        delete_files (pos_files_non)
        delete_files (pos_files_wn)
        
        print "Parsing pileup output to get base counts!\n"
        com_af_files = parse_pileup_output(patient_files, patient_snv_wn, patient_snv_non)
         
        print "Adjusting tumour base count by tumour content!\n"
        for type in com_af_files:
            adjust_allele_frequency(com_af_files[type], patient_files)
            #print com_af_files[type]
        #print com_af_files
        print "Performing Fisher Exact Test for tumour/normal pairs!\n"
        run_fisher_exact_test(Rscript_path, r_script)
        
        
        print "Deleting intermediate files!" 
        #extension=['combined.adjusted', 'combined', '.pileup.log', 'pileup.AFcounts']
        extension=['combined.adjusted', 'combined', '.pileup.log']
        #delete_files_by_extension(extension)
        
        print "Reading in DNA and RNA af files and making af dictionary! \n"
        af_out = make_af_dict(patient_files)
        DNA_af_wn_dict = af_out[0]["DNA"][0]
        DNA_af_non_dict = af_out[0]["DNA"][1]
        RNA_af_wn_dict = af_out[0]["RNA"][0]
        RNA_af_non_dict = af_out[0]["RNA"][1]
        DNA_AFcounts_afs_wn = af_out[1]["DNA"]
        RNA_AFcounts_afs_wn = af_out[1]["RNA"]
        #print "Deleting AFcounts.af files!\n"
        #delete_files(AFcounts_afs_wn)
         
        print "Combining snv_non summary, af, and annotation results!\n"
        pairing_status = 'unpaired'
        combine_sum_af_anno(pairing_status, hm_snv_sum_non, hm_snv_sum_non_tmp, 
                            DNA_af_non_dict, RNA_af_non_dict, patient_files, sum_snv_header_non)
    
        
        print "Combining snv_wn summary, af, and annotation results!\n"
        pairing_status = 'paired'
        combine_sum_af_anno(pairing_status, hm_snv_sum_wn, hm_snv_sum_wn_tmp, 
                            DNA_af_wn_dict, RNA_af_wn_dict, patient_files, sum_snv_header_wn)
        
        print "Deleting intermediate files!"
        #extension=['adjusted.fisher', '.pileup', '.vcf', 'AFcounts.af']
        extension=['adjusted.fisher', '.vcf', 'AFcounts.af']
        #delete_files_by_extension(extension)
        '''


        # filter variants
        variant_summary = hm_snv_sum_wn 
        filtered_summary = filter_variants(variant_summary)
  
        # make expands input files
        make_expands_input_file(filtered_summary, patient_files_wn)
        sys.exit() 
    
    if args['data_type'] == 'spc':
        #out_dict2 = make_patient_vcfs_dict_to_fix(patient_files)
        #patient_wn_snv_vcfs = out_dict2[0]
        #patient_wn_indel_vcfs = out_dict2[1]
        #patient_non_vcfs = out_dict2[2]
        print "mmmmmmmmmmmmm", patient_files 
        '''
        patient_file = "NB_333_patients.tsv"
        vcf_file = "single_mutseq_strelka_vcfs_for_333.tsv"
        patients = make_a_list(patient_file)
        vcfs = make_a_list(vcf_file)
        sum_out_non = summarize_snvs_for_clonality(patients, vcfs, sum_snv_header, hm_snv_sum_non_tmp, low_snv_sum_non_tmp, modifier_snv_sum_non_tmp)
        snv_list = sum_out_non[0]
        snv_pos_file = sum_out_non[1]


        print "Generating mpileup scripts!\n"
        out = make_pileup_scripts(patient_files, snv_pos_file)
        pileup_scripts = out[0]
        pileup_complete_stamps = out[1]
        print pileup_complete_stamps
        #sys.exit()
        
        print "Qsub mpileup scripts!\n"
        #qsub_scripts(pileup_scripts)
 
        print "Detecting if pileup jobs on cluster finised!\n"
        #detect_cluster_jobs(pileup_complete_stamps)
        '''
        '''
        print "Deleting files! \n"
        delete_files(pileup_scripts)
        delete_files(pileup_complete_stamps)
        delete_files (pos_files_non)
        delete_files (pos_files_wn)
        '''
        print "Parsing pileup output to get base counts!\n"
        parse_pileup_output(patient_files, snv_list)
        '''
        print "Deleting intermediate files!"
        extension=['combined.adjusted', 'combined', '.pileup.log', 'pileup.AFcounts']
        delete_files_by_extension(extension)
 
 
        ''' 
        
        print "Reading in af files and making af dictionary! \n"
        af_non_dict = make_af_dict(patient_files)
 
        print "Deleting AFcounts.af files!\n"
        #delete_files(AFcounts_afs_wn)
 
        print "Combining snv_non summary, af, and annotation results!\n"
        pairing_status = 'unpaired'
        combine_sum_af_anno(pairing_status, hm_snv_sum_non, hm_snv_sum_non_tmp,
                            af_non_dict, patient_files, sum_snv_header_non)
        ''' 
        print "Deleting intermediate files!"
        extension=['adjusted.fisher', '.pileup', '.vcf', 'AFcounts.af']
        delete_files_by_extension(extension)
        '''
        print "Making expand input files!\n"
        make_expands_input_file(hm_snv_sum_non, patient_files) 

    print "Summarization scripts finished on: %s\n" % datetime.datetime.now()    
    ##################
    # updated nov 17, 2015, combine paired DNA RNA  results still need to be worked out
    ###########################       
if __name__ == '__main__':
    __main__()

