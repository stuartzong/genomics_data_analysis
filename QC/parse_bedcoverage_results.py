import os
import os.path
import time
import subprocess
import re
import sys
import glob
#import argparse
import csv
from collections import defaultdict
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt


def library_to_path(library):
    return os.path.join('/home/szong/projects/development/coverage', '{0}.csv'.format(library))


def get_data(libraries):
    """Read bedcoverage data for each position from CSV files."""
    df = pd.DataFrame()
    for library in libraries:
        df_temp = pd.read_csv(library_to_path(library), header=None, index_col=[4],
                              usecols=[5], na_values=['nan'])
        # df_temp = df_temp.rename(columns={: library})
        df = df.join(df_temp)
    return df    


def library_list(infile):
    libraries = []
    with open(infile, 'r') as fh:
        for line in fh:
           libraries.append(line.strip()) 
    print libraries
    return libraries


def make_numreads_dict(bamstats_summary):
    num_reads_dict = dict()
    with open(bamstats_summary) as fh:
        for line in fh:
            print line
            sl = line.rstrip().split('\t')
            lib = sl[0]
            num_reads = sl[1]
            num_reads_dict[lib] = num_reads
        print num_reads_dict
        return num_reads_dict
  

def make_exon_dict(exon_file):
  with open (exon_file) as fh:
      exon_d = {}
      exon_len_d = {}
      for line in fh:
         sl = line.rstrip().split('\t')
         start = sl[1]
         end = sl[2]
         exon_len = int(end) - int(start)
         exon  = "_".join(["exon", start, end])
         exon_len_d[exon] = exon_len
         try:
           exon_d[exon] = [start, end]
         except KeyError:
           print "KeyError!"
      print exon_d
      print exon_len_d
      return [exon_d, exon_len_d]



def normalize_coverage(reads_dict, exons, normalized_file, coverage_file):
    #with open(normalized_file,  'wb') as writer:
    #writer = csv.writer(opf,  delimiter='\t')
    # write as per exon
    for exon in exons:
        exon_file = "".join(["exon_", str(exon), "_coverage.txt"])
        exon_start = int(exons[exon][0])
        exon_end = int(exons[exon][1])
        with open(exon_file,  'wb') as writer2:
            with open(coverage_file,  'r') as fh:
                records = csv.DictReader(fh,  delimiter='\t')
                headers = records.fieldnames
                libraries = []
                print headers
                for header in sorted(set(headers)):
                    if ("A" in header):
                         print header
                         libraries.append(header)
                header_start = ["chr", "start", "end", "exon", "pos"]	
                #writer.write("\t".join(header_start + libraries))
                #writer.write("\n")
                #make a list of aligned reads
                aligned_reads = []
                for library in libraries:
                    aligned_reads.append(int(reads_dict[library]))
                average_reads = sum(aligned_reads)/len(aligned_reads)
                print "xxxxxxx", average_reads 
                writer2.write("\t".join(header_start + libraries))
                writer2.write("\n")
                for line in records:
                    start = line["start"]
                    pos = line["pos"]
                    coord = int(start) + int(pos) 
                    #print "bbbb", start, pos, coord, exon_start, exon_end
                    if (coord >= exon_start and coord <= exon_end):
                        coverages = []
                        for library in libraries:
                            coverages.append(line[library])
                        normalized_cov = [str(int(average_reads*float(coverage) / int(aligned_read))) for coverage, aligned_read  in zip(coverages, aligned_reads)]
                        print normalized_cov
                        positions = []
                        for i in header_start:
                            positions.append(line[i])
                        writer2.write("\t".join(positions + normalized_cov))
                        writer2.write("\n")
           

def normalize_coverage_per_base(reads_dict, coverage_file):
    normalized_file = ".".join([coverage_file, "normalized"])
    with open(normalized_file,  'wb') as opf:
            writer = csv.writer(opf,  delimiter='\t')
            with open(coverage_file,  'r') as fh:
                records = csv.DictReader(fh,  delimiter='\t')
                headers = records.fieldnames
                libraries = headers[5:]
                header_start = ["chr", "start", "end", "exon", "pos"]   
                print libraries
                aligned_reads = [float(reads_dict[i]) for i in libraries]
                average_reads = numpy.mean(aligned_reads)
                writer.writerow(headers)
                for line in records:
                    leading_columns = [line[i] for i in header_start]
                    coverages = [float(line[lib]) for lib in libraries]
                    print aligned_reads
                    print coverages
                    normalized_cov = [str(float(average_reads)*float(coverage) / float(aligned_read)) for coverage, aligned_read  in zip(coverages, aligned_reads)]
                    content = leading_columns + normalized_cov
                    writer.writerow(content)
                print libraries, aligned_reads, average_reads, coverages 
    return average_reads

def make_lib_pat_dict(lib_pat_file):
    #patlib_d patient_library_dictionary
    patlib_d = dict()
    with open (lib_pat_file) as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            lib = line['lib']
            pat = "_".join([line['patient'].split("-")[2], line['status']])
            patlib_d[pat] = lib
        print patlib_d
        return patlib_d


def identify_mutants(exons, patlib_d, summary_file):
    for exon in exons:
        exon_file = "".join(["exon_", str(exon), "_mutants.txt"])
        exon_start = int(exons[exon][0])
        exon_end = int(exons[exon][1])
        with open(exon_file,  'wb') as writer:
            with open(summary_file,  'r') as fh:
                records = csv.DictReader(fh,  delimiter='\t')
                for line in records:
                    gene = line['Gene']
                    if ( "TP53" in gene):
                        chr = line['Chromosome']
                        pos = int(line['position'])
                        pat = line['PatientID']
                        if (pos >= exon_start and pos <= exon_end):
                            lib = patlib_d[pat]
                            writer.write(lib)
                            writer.write("\n")
def make_lib_pathology_names_dict(infile):
    lib_pathology_d = dict()
    with open (infile) as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        for line in records:
            lib = line['library']
            pat = line['patient']
            pathology = line['pathology']
            value = ":".join([pat, pathology])
            lib_pathology_d[lib] = value
        print lib_pathology_d
        return lib_pathology_d

def add_info(reads_dict, lib_pathology_d, exon_cov_file, ave_aligned_reads, exon_len_d):
    out_file = ".".join([exon_cov_file, "final"])
    with open(out_file,  'wb') as opf:
            writer = csv.writer(opf,  delimiter='\t')
            with open(exon_cov_file,  'r') as fh:
                records = csv.DictReader(fh,  delimiter='\t')
                headers = records.fieldnames
                extra_header = ["patient", "pathology", "aligned_reads"] + [i for i in headers if (("exon" in i) or ("deviation" in i))]
                writer.writerow(headers + extra_header)
                for line in records:
                    lib = line['lib_name']

                    patient = lib_pathology_d[lib].split(':')[0]
                    pathology = lib_pathology_d[lib].split(':')[1]
                    aligned_reads = reads_dict[lib]
                    line_content = [line[i] for i in headers]
                    adj_factor = float(ave_aligned_reads)/ int(aligned_reads)
                    
                    extra_content = [float(line[i].rstrip())*adj_factor for i in headers if (("exon" in i) or ("deviation" in i))]
                    exon_lens = [exon_len_d[i] for i in headers if ("exon"  in i)]
                    extra_exon_cov = [float(line[i].rstrip())*adj_factor for i in headers if ("exon" in i)]
                    extra_std_avg = numpy.mean([float(line[i].rstrip())*adj_factor for i in headers if ("deviation" in i)])
                    
                    total_exon_len = numpy.sum(exon_lens)
                    all_exon_total = numpy.sum([int(exon_lens[i])*float(extra_exon_cov[i]) for i in range(len(exon_lens))])
                    all_exon_avg = all_exon_total/total_exon_len
                    content = line_content + [patient, pathology, aligned_reads] + extra_content + [all_exon_avg, extra_std_avg]
                    writer.writerow(content)
    
 
def __main__():
    # wkdir = os.getcwd()
    # bamstats_summary = "bamstats_summary.txt"
    library_list = 'library_list.txt'
    libraries = library_list(library_list)
    get_data(libraries)
    sys.exit()
    #original coverage file
    coverage_file = "combined_exon_coverage.txt"
    exon_file = "WT1_intervals.txt"
    exon_len_d = make_exon_dict(exon_file)[1]
    # sys.exit()
    reads_dict = make_numreads_dict(bamstats_summary)
    ave_aligned_reads = normalize_coverage_per_base(reads_dict, coverage_file)
    
    sys.exit()
    lib_pathology_file = "WT_pathology.csv"
    lib_pathology_d = make_lib_pathology_names_dict(lib_pathology_file)
    exon_cov_file = "WT_MYCN_exon_coverage_and_length_of_intron_covered.csv" 
    add_info(reads_dict, lib_pathology_d, exon_cov_file, ave_aligned_reads, exon_len_d)


    sys.exit()
    lib_pat_file = "vcf_bam_path.txt"
    patlib_d = make_lib_pat_dict(lib_pat_file)
    
    summary_file = "SNV_summary_no_normal.txt.filtered.final"
    identify_mutants(exons, patlib_d, summary_file)
if __name__ == '__main__':
    __main__()
