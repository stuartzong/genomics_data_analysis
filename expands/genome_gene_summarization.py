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
1. Prepare vcf files for Bedtools Bedintersect
2. Overlap regions in vcfs with specific targeted regions
3. Remove all off-target variants called for all patients
4. Summarize somatic SNVs and INDELS events called in vcf at gene and variant level
5. Run samtools fileup for all the SNV positions
6. Parse pileup output to get the base count for each SNV position
7. Perform Fisher Exact Test to calculate p value for base changes
8. Annotate SNVs with dbSNPS, COSMIC, Ensembl full gene name, and Strelka high confidence calls

Warnings:
Hardcoded paths present in this version.

Input files:
1. Tumour content file
2. Paths for raw vcf files
3. Ensembl gene name file
4. Cosmic annotation file
5. Paths to Strelka vcf files 

"""

def make_patient_vcf_file_dict(vcf_file):
    #initiate a dictionary, patient:status:[lib, vcf_file, bam_file, strelka_vcf, tc]
    patient_status_files_dict = dict()
    vcfs = []
    #print "The file is: %s.\n" % vcf_file
    with open (vcf_file, 'r') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        tempDict=dict()
        for line in records:
            patient = line['patient']
            status = line['status']
            maf_path = line['maf_path']
            value = [maf_path]
            vcfs.append( maf_path )
            try:
                patient_status_files_dict[patient][status] = value
            except KeyError:
                if patient not in patient_status_files_dict:
                    patient_status_files_dict[patient] = {}
                if status not in patient_status_files_dict[patient]:
                    patient_status_files_dict[patient][status] = value
            
        return [patient_status_files_dict, vcfs ]
  


      
def run_intersectBed(fileA, fileB, patient):
    #run bedtools intersectBed
    """use subprocess to run shell command"""
    bedtoolsPath = "/home/rcorbett/aligners/bedtools/BEDTools-Version-2.15.0/bin"
    intersectedFile = ".".join([patient, "intersected.vcf"])
    p = subprocess.Popen('%s/intersectBed -a %s -b %s -wa > %s' % (bedtoolsPath, fileA, fileB, intersectedFile),  shell=True, stdout=subprocess.PIPE)
    output,  err = p.communicate()
    #print output
   

def list_files_by_extension(extension, vcf_fileList):
    filelist = [ f for f in os.listdir(".") if f.endswith(tuple(extension))]
    writer = open(vcf_fileList, 'w')
    for f in filelist:
        #print "writing vcf files into vcf_fileList: %s" % f
        writer.write(f)
        writer.write("\n")


def delete_files_by_extension(extension):
    #clean up directory
    filelist = [ f for f in os.listdir(".") if f.endswith(tuple(extension))]
    for f in filelist:
        print "Removing files: %s" % f
        os.remove(f)


def make_patient_status_files_nested_dict(infile):
    """
    read in a file,  output a dictionary of dictionay,  in this case,  read in\
    the patient-primary/relapse lib names, vcf_path, and bam path.
    """
    newDict=dict()
    with open(infile,  'r') as handle:
        records = csv.DictReader(handle,  delimiter='\t')
        tempDict=dict()
        for line in records:
            key1 = line['patient']
            key2=line['status']
            value=[line['lib'], line['vcf_path'], line['bam_path']]
            tempDict[key2]=value
            #print "reading:", key1, key2, value
            try:
                newDict[key1][key2]=value
            except KeyError:
                newDict[key1]=dict(tempDict)
        #print newDict
        return newDict
  
  


def summarize_snvs(vcf_files, summaryFile):
    tree = lambda: defaultdict(tree)
    geneSnvPatientDict = tree()
    geneDict = dict()
    snvDict = dict()
    with open(vcf_files, 'r') as fh:
       for f in records:
          vcf_file = f.strip()
          with open(vcf_file,  'r') as handle:
              records = csv.DictReader(handle,  delimiter='\t')
              for line in records:
                  chr, start, end = sl[0], sl[1], sl[2]
                  chr = line['']
                  start = line['']
                  end = line['']
                  # CGI maf has overlapped gene names
                  multi_genes = sl[3].split('|')
                  gene_list = [gene for gene in multi_genes if gene in targetted_genes]              
                  try:
                      gene = gene_list[0]
                  except:
                      print "<<<<<<<<<<<<<<<<<"
                      print multi_genes
                      print gene_list
                      gene = multi_genes[0]
                      print gene                  
                  impact_type, cosm_id = sl[6], sl[14]
                  ref_base, tum_allele1, tum_allele2, nor_allele1, nor_allele2 = sl[16], sl[17], sl[18], sl[19], sl[20]
                  tum_alt_count, tum_ref_count, tum_total = sl[8], sl[9], sl[10]
                  nor_alt_count, nor_ref_cout, nor_toatl = sl[11], sl[12], sl[13]

                  snv = ":".join([chr, start, ref_base, tum_allele2])
                  patient = re.split('[-.]', vcf_file)[2]
                  dbSNPid = sl[7]
                  COSMICid = sl[14]
 
                  #make gene>patient dictionary to count how many patients have variant in this gene
                  try:
                     geneDict[gene].append(patient)
                  except:
                     geneDict[gene] = [patient]
                  #get gmaf value      
                  details = "_".join([dbSNPid, COSMICid, impact_type, tum_alt_count, tum_ref_count, tum_total, nor_alt_count, nor_ref_cout, nor_toatl])
                  
                  #put all relevant info into a triple nested dictionary: gene>snv>patient>details
                  try:
                     geneSnvPatientDict[gene][snv][patient].append(details)
                  except:
                     geneSnvPatientDict[gene][snv][patient]=[details]
       #print "******************geneSnvPatientDict is:", geneSnvPatientDict
       #write info into summary file
       writer = open(summaryFile, 'w')
       writer.write("\t".join(["Gene",  "numPatientGeneLevel", "PatientIDsGeneLevel", "numSNVsForThisGene", "Chromosome", "position", "referenceBase", "AlternativeBase", "numPatientSnvLevel", "PatientIDsSnvLevel", "PatientID", "dbSNPid",  "COSMICid", "impactDetails", "tumourAlternativeBaseCount", "tumourReferenceBaseCount", "tumourSequencingCoverage", "normalAlternativeBaseCount", "normalReferenceBaseCount", "normalSequencingCoverage\n"]))
       for gene in geneSnvPatientDict:
         numPatientGeneLevel = len(list(set(geneDict[gene])))
         numSNVs = len (geneSnvPatientDict[gene].keys())
         for snv in geneSnvPatientDict[gene]:
           numPatientSnvLevel = len(geneSnvPatientDict[gene][snv])
           snvPatient = []
           snvDetails = []
           for patient in geneSnvPatientDict[gene][snv]:
              snvDetails.append(details)
              snvPatient.append(patient)
           for patient in geneSnvPatientDict[gene][snv]:
              geneLevelPatient = ",".join(list(set(geneDict[gene])))      
              svnLevelPatient = ",".join(list(set(snvPatient)))        
              #snvDetails: 10:64927823:C:G
              snvDetails = snv.split(":")              
              #print "XXXXXXXX",geneSnvPatientDict[gene][snv][patient][0]
              snvAnnoDetails = list(set(geneSnvPatientDict[gene][snv][patient]))[0].split("_")
              writer.write("\t".join([gene,  str(numPatientGeneLevel), geneLevelPatient, str(numSNVs), snvDetails[0], snvDetails[1], snvDetails[2], snvDetails[3], str(numPatientSnvLevel), svnLevelPatient, patient, snvAnnoDetails[0], snvAnnoDetails[1], snvAnnoDetails[2],snvAnnoDetails[3],snvAnnoDetails[4],snvAnnoDetails[5], snvAnnoDetails[6],snvAnnoDetails[7], snvAnnoDetails[8]]))
              writer.write("\n")

    return geneSnvPatientDict

def summarize_indels(vcf_files, summaryFile):
    tree = lambda: defaultdict(tree)
    geneIndelPatientDict = tree()
    geneDict = dict()
    indelDict = dict()
    with open(vcf_files, 'r') as fh:
       for f in fh:
          vcf_file = f.strip()
          with open(vcf_file,  'r') as handle:
              for line in handle:
                  sl = line.strip().split("\t")
                  Info = sl[7]
                  #print "Info is:",  Info
                  #only look for high and moderate impact event,  excluding INDELs
                  if (("HIGH" in Info) or ("MODERATE" in Info)) and (Info.startswith("INDEL")):
                     splitInfo = Info.split(";")
                     transcripts =  [i for i in splitInfo if i.startswith("EFF=")][0].split("=")[1].split(",")
                     HM_ImpactTranscripts = [k for k in transcripts if ("HIGH" in k or "MODERATE" in k)]
                     #pick the first transcript to get the gene name
                     gene = HM_ImpactTranscripts[0].split("|")[5]
                     transcriptDetails = HM_ImpactTranscripts[0]
                     indel = ":".join([sl[0], sl[1], sl[3], sl[4]])
                     patient = re.split('[-.]', vcf_file)[2]




                     dbSNPid_tmp = sl[2].split(";")[0]
                     if dbSNPid_tmp.startswith("rs"):
                        dbSNPid = dbSNPid_tmp
                     else:
                        dbSNPid = "novelSNP"
                     
                     try:
                        COSMICid = sl[2].split(";")[1]
                     except:
                        if dbSNPid_tmp.startswith("COS"):
                            COSMICid = dbSNPid_tmp
                        else:
                            COSMICid = "not_in_COSMIC"


                     #make gene>patient dictionary to count how many patients have variant in this gene
                     try:
                        geneDict[gene].append(patient)
                     except:
                        geneDict[gene] = [patient]
                     #get gmaf value      
                     try:
                        gmaf = [j for j in splitInfo if j.startswith("GMAF=")][0].split("=")[1]
                     except:
                        gmaf = "gmaf_unknown"
                     details = ":".join([dbSNPid, COSMICid, gmaf, transcriptDetails])

                     #put all relevant info into a triple nested dictionary: gene>indel>patient>details
                     try:
                        geneIndelPatientDict[gene][indel][patient].append(details)
                     except:
                        geneIndelPatientDict[gene][indel][patient]=[details]
       #print "******************geneIndelPatientDict is:", geneIndelPatientDict
       #write info into summary file
       writer = open(summaryFile, 'w')
       writer.write("\t".join(["Gene",  "numPatientGeneLevel", "PatientIDsGeneLevel", "numINDELsForThisGene", "Chromosome", "position", "referenceBase", "AlternativeBase", "numPatientIndelLevel", "PatientIDsIndelLevel", "PatientID", "dbSNPid", "GMAF", "COSMICid", "ImpactDetails\n"]))
       for gene in geneIndelPatientDict:
         numPatientGeneLevel = len(list(set(geneDict[gene])))
         numINDELs = len (geneIndelPatientDict[gene].keys())
         for indel in geneIndelPatientDict[gene]:
           numPatientIndelLevel = len(geneIndelPatientDict[gene][indel])
           indelPatient = []
           indelDetails = []
           for patient in geneIndelPatientDict[gene][indel]:
              indelDetails.append(details)
              indelPatient.append(patient)
           for patient in geneIndelPatientDict[gene][indel]:
              geneLevelPatient = ",".join(list(set(geneDict[gene])))      
              svnLevelPatient = ",".join(list(set(indelPatient)))        
              #indelDetails: 10:64927823:C:G
              indelDetails = indel.split(":")              
              indelAnnoDetails = list(set(geneIndelPatientDict[gene][indel][patient]))[0].split(":")
              writer.write("\t".join([gene,  str(numPatientGeneLevel), geneLevelPatient, str(numINDELs), indelDetails[0], indelDetails[1], indelDetails[2], indelDetails[3], str(numPatientIndelLevel), svnLevelPatient, patient, indelAnnoDetails[0], indelAnnoDetails[2], indelAnnoDetails[1], indelAnnoDetails[3]]))
              writer.write("\n")

    return geneIndelPatientDict




def pre_filter_maf(patient, vcf_file, maf_out_file):
  #patient = vcf_file.strip().split("/")[8]
  #maf_out_file = ".".join([patient, "maf.filtered"])
  with open(vcf_file, 'rU') as handle:
       print "processing ", vcf_file
       records = csv.DictReader(handle,  delimiter='\t')
       #outFile = "{patient}.tmp".format(patient=patient)
       with open(maf_out_file, 'wb') as writer:
         for line in records:
             try:
                #print line
                var_type  = line['VariantType']
                if (var_type == ""):
                    var_type = "NA"
             
                mut_status = line["Mutation_Status"]
                if (mut_status == ""):
                    mut_status = "NA"
                chr = line['Chromosome']
                start = line['Start_position']
                end = line['End_position']
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
                #t_af = float(tum_alt)/float(tum_total)
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
                nor_allele1 = line['Match_Norm_Seq_Allele1']	
                nor_allele2 = line['Match_Norm_Seq_Allele2']	
                if (nor_allele2 == ""):
                    nor_allele2 = "N"
                tum_barcode = line['Tumor_Sample_Barcode']
                #print mut_status, var_type
                #if (mut_status == "Somatic") and ((var_type == "SNP") or (var_type == "Del") or (var_type == "Ins") ):
                if (mut_status == "Somatic") and (var_type == "SNP"): # and (int(tum_total) >= 10) and  (int(nor_total) >= 10) and ( int(nor_alt) < 1) and (t_af > 0.05):
                    writer.write("\t".join([chr, start, end, gene, var_type, mut_status,
                                           var_classification, dbsnp, tum_alt, tum_ref, tum_total,
                                           nor_alt, nor_ref, nor_total, cosm_id, cosm_gene,
                                           ref, tum_allele1, tum_allele2, nor_allele1, nor_allele2]))
                    writer.write("\n") 
         
             except:#there are apparent mistakes in the CGI maf files
                 print "Invalid entry in original CGI maf file skipped!\n", 
                 print line
                 continue
    
    
    
    

def __main__():
    input_files = "CGI_292_vcfs.csv"
    #patient:status:[lib, vcf_file, bam_file, strelka_vcf]
    out = make_patient_vcf_file_dict(input_files)
    patient_files_dict = out[0]
    vcfs = out[1]

    #remove vcf header and run intersectBed to remove off target regions
    for patient in patient_files_dict:
       for status in patient_files_dict[patient]:
           ''' 
           try:
               vcf_file = patient_files_dict[patient]["Primary"][0]
           except KeyError:
               vcf_file = patient_files_dict[patient]["Relapse"][0]
           '''   
           vcf_file = patient_files_dict[patient][status][0]
           maf_out_file = ".".join([patient, status])
           pre_filter_maf(patient, vcf_file, maf_out_file)      

    '''
    #make intersected vcf file list:
    extension = ['.intersected.vcf']
    vcf_fileList = "vcf.files"
    list_files_by_extension(extension, vcf_fileList)

    vcf_files = vcf_fileList

    pprint( patient_files_dict)
    snvSummaryFile = "SNV.summary.txt"
    summarize_snvs(vcfs, snvSummaryFile)
    #indelSummaryFile = "INDEL.summary.txt"
    #summarize_indels(vcf_files, indelSummaryFile)
    '''
if __name__ == '__main__':
    __main__()

