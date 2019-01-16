import collections
import headers as HEADER
import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from pprint import pprint
from itertools import islice
import operator
import ConfigParser
from collections import defaultdict
import headers as HEADER
import numpy


from jinja2 import Environment, FileSystemLoader
import logging
import colorlog

logger = colorlog.getLogger()
logger.setLevel(logging.DEBUG)
handler = logging.StreamHandler()
handler.setLevel(logging.DEBUG)
handler.setFormatter(
    colorlog.ColoredFormatter('%(log_color)s%(levelname)s:%(name)s:%(message)s'))
logger.addHandler(handler)


"""
To run this script, type the following command:
/gsc/software/linux-x86_64/python-2.7.2/bin/python summarize_variants_newest_version.py -dt spc > summary_log_file.txt

"""

def get_files(bam_vcf_files, input_headers):
    """ 
    Dictionary: holding all file paths
    {patient ->
             {status ->
                     {file_identifier -> file_path}}}  
    """
    patient_files = dict()
    with open(bam_vcf_files, 'r') as fh:
        records = csv.DictReader(fh,  delimiter='\t')
        headers = set(records.fieldnames)
        input_headers = set(input_headers)
        # check if all mandatory columns prsent in input file
        if input_headers.issubset(headers):
            logger.info("input file has all mandatory columns! Continue...")
            for line in records:
                patient = line['patient']
                status = line['status']
                if patient not in patient_files:
                    patient_files[patient] = {}
                if status not in patient_files[patient]:
                    patient_files[patient][status] = line
                else:
                    logger.error('Duplicate status entries!'\
                                 'One tissue sequenced multiple times?')
                    sys.exit()
        else:
            logger.critical("Input file doesn't contain all mandatory headers!\n %s"
                            % list(input_headers))
            sys.exit()
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
        exist_status = False
        read_status = False
        exist_status = exist(file)
        if (exist_status):
            read_status = is_other_readable(file)
        if (not (exist_status and read_status)):
            missing_files.append(file)
    if not missing_files:
        print "All files exist and are readable!\n" 
    else:
        print "ERROR: The following files are either missing or no read permission!\n"
        print "ATTENTION: Please check to make sure STRELKA has completed for this sample!\n"
        for file in missing_files:
            print file
        ##sys.exit()


def make_BEDcoverage_script(patient_files, bed_regions, template_dir):
    bed_path = "/home/rcorbett/aligners/bedtools/BEDTools-Version-2.15.0/bin/"
    sam_path = "/home/rcorbett/aligners/samtools/samtools-0.1.17/"
    wkdir = os.getcwd()
    scripts = []
    stamps = []
    coverage_files = []
    for patient in patient_files:
        for status in patient_files[patient]:
            library = patient_files[patient][status]['DNA_lib']
            bam = patient_files[patient][status]['DNA_bam'] 
            print "xxxxxxxxxxxx", bam
            complete_stamp = ".".join([patient, status, library, "complete"])
            stamps.append(complete_stamp)
            coverage_file = ".".join([patient, status, library, "exon.coverage.txt"])
            coverage_files.append(coverage_file)
            script = populate_template(patient,
                                       status,
                                       wkdir,
                                       library,
                                       bed_regions,
                                       bam,
                                       template_dir)
 
            scripts.append(script)
    return [scripts, stamps, coverage_files]
    
    
def populate_template(patient, status, wkdir, library,
                      bed_regions, bam, template_dir):
    bedcoverage_script = ".".join([patient, status, library, "BEDcoverage.sh"])
    jinja2_env = Environment(loader=FileSystemLoader([template_dir]),
                             trim_blocks=True)
    template = jinja2_env.get_template('bedcoverage_template.sh')
    with open(bedcoverage_script, 'wb') as opf:
        content = template.render(patient=patient,
                                  status=status,
                                  library=library,
                                  bam=bam,
                                  bed_regions=bed_regions,
                                  wkdir=wkdir)
        opf.write(content)
        logger.info('templated {0}'.format(bedcoverage_script))
    return bedcoverage_script


def delete_files_by_extension(extension):
    files = [ f for f in os.listdir(".") if f.endswith(tuple(extension))]
    for f in files:
        os.remove(f)

def delete_files(files):
    for f in files:
        os.remove(f)

  
def make_lib_bam_dict(infile):
    # {'A44623': ['/path/aa.bam'], ...} 
    d = dict()
    with open(infile, 'r') as fh:
        for line in fh:
           sl = line.rstrip().split('\t')
           lib = sl[0]
           bam = sl[1]
           try:
             d[lib].append(bam)
           except KeyError:
             d[lib] = [bam]
        return d

def qsub_scripts(scripts):
    """ qsub scripts """
    wkdir = os.getcwd()
    for script in scripts:
        p = subprocess.Popen('ssh tachpc \"cd %s;  qsub %s\"' % (wkdir, script),  shell=True, stdout=subprocess.PIPE)
        output,  err = p.communicate()


def detect_cluster_jobs(complete_stamps):
    """ detect if job on cluster finised """
    completed = False
    print "Waiting for cluster jobs to finish!\n"
    while (not completed):
        time.sleep(10)
        for file in complete_stamps:
            if (os.path.exists(file)):
                completed = True
            else:
                completed = False
                break
    print "All cluster jobs finished? %s" % completed
 

def combine_coverage(coverage_files, combined_file, header):
    d = dict()
    for coverage_file in coverage_files:
        lib = '_'.join(coverage_file.split(".")[:3])
        with open(coverage_file, 'r') as fh:
            for line in fh:
                sl = line.split("\t")
                key = "_".join(sl[:5])
                cov = sl[5]
                if lib not in d:
                    d[lib] = {}
                if key not in d[lib]:
                    d[lib][key] = cov 
    libs = [lib for lib in sorted(d)]
    header_final = header + libs
    print "bbb", header_final

    with open(combined_file, 'wb') as opf:
        writer = csv.writer(opf,  delimiter='\t')
        writer.writerow(header_final)
        # access the first key in d
        first_lib = next (iter (d.keys()))
        print "aaaaaaa", first_lib
        keys = [key for key in d[first_lib]]
        #lambda is an anonymous function, sort string_num mix on numerical order
        keys = sorted(keys, key=lambda item: (int(item.split('_')[4])
                               if item[4].isdigit() else float('inf'), item))
        for key in keys:
            covs = [d[lib][key].replace('\n', '') for lib in libs]
            content = key.split('_') + covs
            writer.writerow(content)

def make_exon_intron_dict(exon_file):
    # output orderedDict exon_d, and intron_d
    with open(exon_file, 'r') as fh:
        positions = []
        exon_d = dict()
        intron_d = dict()
        for line in fh:
            sl = line.split('\t')
            start = sl[1]
            end = sl[2]
            positions.append(start)
            positions.append(end)
            key = "_".join(["exon", start, end])
            exon_d[key] = [start, end]
        first = int(positions[0]) - 1000
        last = int(positions[-1]) + 1000
        positions.append(str(last))
        positions = [first] + positions
        #print positions
    # make intron dict
    for i in range(0,len(positions)/2):
        start_idx = 2*i
        end_idx = start_idx + 1
        #print start_idx, end_idx
        start = str(positions[start_idx])
        end = str(positions[end_idx])
        key = "_".join(["intron", start, end])
        intron_d[key] = [start, end]
    exon_d = collections.OrderedDict(sorted(exon_d.items()))
    intron_d = collections.OrderedDict(sorted(intron_d.items()))
    print "intron dict is: \n"
    pprint(intron_d)
    print "exon dict is:\n"
    pprint(exon_d)
    return [exon_d, intron_d]

def calculate_exon_coverage_to_delete(BEDcoverage_file, exon_dict, intron_dict, cov_dict_all, first_exon_pos):
    start_offset = first_exon_pos - 1000 
    cov_dict = dict()
    ontron_cov_dict = dict()
    ave_cov_dict = dict()
    print "Proccessing %s\n" % BEDcoverage_file
    with open(BEDcoverage_file, 'r') as fh:
        for line in fh:
            line = line.rstrip('\n')
            sl = line.split('\t')
            position = int(sl[-2]) + start_offset
            coverage = int(sl[-1])
            # bin coverages as per exon
            for exon in exon_dict:
                start = int(exon_dict[exon][0])
                end = int(exon_dict[exon][1])
                if (position in range(start, end)):
                    try:
                        cov_dict[exon].append(coverage)
                    except KeyError:
                        cov_dict[exon] = [coverage]
                    try:
                        cov_dict_all[exon].append(coverage)
                    except KeyError:
                        cov_dict_all[exon] = [coverage]

           # bin coverages as per intron
            for intron in intron_dict:
                start = int(intron_dict[intron][0])
                end = int(intron_dict[intron][1])
                if (position in range(start, end)):
                    try:
                        intron_cov_dict[intron].append(coverage)
                    except KeyError:
                        intron_cov_dict[intron] = [coverage]
                    
    
    ave_cov_dict = calculate_mean_std_covered_intron(cov_dict, intron_cov_dict)
    print intron_cov_dict
    return [ave_cov_dict, cov_dict, cov_dict_all]
                 
def calculate_mean_std_covered_intron(cov_dict, intron_cov_dict):
    ave_exintron_cov_dict = dict()
    for exon in cov_dict:
        # average + std
        cov_list = cov_dict[exon]
        average = numpy.mean(cov_list)
        std = numpy.std(cov_list)
         
        # length of covered intron
        left_intron = exon.split('_')[1]
        right_intron = exon.split('_')[2]
        print "yyyyyyyy", exon, left_intron, right_intron
        print "intron_cov_dict is:", intron_cov_dict
        for intron_key in intron_cov_dict:
            if (left_intron in intron_key):
                left_intron_cov = intron_cov_dict[intron_key]
                left_intron_len = len(left_intron_cov)
                left_cov_status = []
                print "llllllleft", left_intron_cov
                for j in range(5):
                    if (int(left_intron_cov[-j])< 10):
                        low_cov = True
                    else:
                        low_cov = False
                    print 
                    left_cov_status.append(low_cov)
                if (False in left_cov_status):
                    for i in range(5, len(left_intron_cov)-5):
                        if (int(left_intron_cov[-i])< 10):
                            low_cov = True
                        else:
                            low_cov = False
                        left_cov_status.append(low_cov)
                        left_cov_status = left_cov_status[1:]
                        if (False not in left_cov_status):
                            left_covered_len = i + 1
                            break
                        else:
                            left_covered_len = left_intron_len
                else:
                    left_covered_len = 0       
                print "left covered len is:", exon, intron_key, left_covered_len, left_intron_len         

            if (right_intron in intron_key):
                right_intron_cov = intron_cov_dict[intron_key]
                right_intron_len = len(right_intron_cov)
                right_cov_status = []
                
                print "rrrrrrrright intron coverage", intron_key, right_intron_cov
                # if 5 continous bases have a coverage less than 10, we deem coverage ends
                for j in range(5):
                    if (int(right_intron_cov[j])< 10):
                        low_cov = True
                    else:
                        low_cov = False
                    print 
                    right_cov_status.append(low_cov)
                if (False in right_cov_status):
                    for i in range(5, len(right_intron_cov)-5):
                        if (int(right_intron_cov[i])< 10):
                            low_cov = True
                        else:
                            low_cov = False
                        right_cov_status.append(low_cov)
                        right_cov_status = right_cov_status[1:]
                        if (False not in right_cov_status):
                            right_covered_len = i + 1
                            break
                        else:
                            right_covered_len = right_intron_len
                else:
                    right_covered_len = 0       
                print "right covered len is:", exon, intron_key, right_covered_len, right_intron_len         
               
        print "kkkkkkkk", exon, left_covered_len, left_intron_len, right_covered_len, right_intron_len
        ave_exintron_cov_dict[exon] = [left_covered_len, left_intron_len, average, std, right_covered_len, right_intron_len]
    return ave_exintron_cov_dict


def calculate_mean_std(cov_dict):
    ave_cov_dict = dict()
    for exon in cov_dict:
        # average + std
        cov_list = cov_dict[exon]
        average = numpy.mean(cov_list)
        std = numpy.std(cov_list)
        ave_cov_dict[exon] = [average, std]
    return ave_cov_dict


#def calculate_covered_introns(exon_dict, intron_dict):
    
 
def calculate_avecov_all_libs(coverage_files, exon_dict, intron_dict, first_exon_pos):
    #calculate per exon coverage
    ave_cov_oufile = "ave_cov_outfile.txt"
    with open (ave_cov_oufile, 'wb') as opf:
        writer = csv.writer(opf,  delimiter='\t')
        cov_dict_all = dict()

        ###############
        start_offset = first_exon_pos - 1000
        # infor holds left, len, exon_ave, righ, len for all exons
        infos = []
        all_lib_ave = []
        for BEDcoverage_file in coverage_files:
            #BEDcoverage_file = "A44754.exon.coverage.txt"
            cov_dict = dict()
            intron_cov_dict = dict()
            ave_cov_dict = dict()

            lib = BEDcoverage_file.split('.')[0]
            print "Proccessing %s\n" % BEDcoverage_file
            with open(BEDcoverage_file, 'r') as fh:
                for line in fh:
                    line = line.rstrip('\n')
                    sl = line.split('\t')
                    position = int(sl[-2]) + start_offset
                    coverage = int(sl[-1])
                    # bin coverages as per exon
                    for exon in exon_dict:
                        start = int(exon_dict[exon][0])
                        end = int(exon_dict[exon][1])
                        if (position in range(start, end)):
                            try:
                                cov_dict[exon].append(coverage)
                            except KeyError:
                                cov_dict[exon] = [coverage]
                            try:
                                cov_dict_all[exon].append(coverage)
                            except KeyError:
                                cov_dict_all[exon] = [coverage]
        
                    # bin coverages as per intron
                    for intron in intron_dict:
                        start = int(intron_dict[intron][0])
                        end = int(intron_dict[intron][1])
                        # print intron, start, end, position, int(sl[-2])
                        if (position in range(start, end)):
                            # print "Yes"
                            try:
                                intron_cov_dict[intron].append(coverage)
                            except KeyError:
                                intron_cov_dict[intron] = [coverage]
            # print "444444444444444444444444", BEDcoverage_file
            # sys.exit()
            # pprint(intron_cov_dict)
            # print "444444444444444444444444", BEDcoverage_file
        
        
            ave_exintron_cov_dict = calculate_mean_std_covered_intron(cov_dict, intron_cov_dict)
            #print "ppppppppp", ave_exintron_cov_dict

            cov_std_covd_intron = [ave_exintron_cov_dict[exon] for exon in sorted(ave_exintron_cov_dict)]
            aa = "_".join(["_".join([str(j) for j in i]) for i in cov_std_covd_intron])
            aa_sl = aa.split('_')
            infos.append(aa_sl)
            exons = [["left_intron_covered", "left_intron_len", exon, "standard_deviation", "right_intron_covered", "right_intron_len"] for exon in sorted(ave_exintron_cov_dict)]
            content = [lib] + aa_sl
            writer.writerow(content)
        exons_aa = ":".join([":".join([str(j) for j in i]) for i in exons])
        ave_cov_all_libs = calculate_mean_std(cov_dict_all)
        cov_std = [ave_cov_all_libs[exon] for exon in sorted(ave_cov_all_libs)]
        bb = "_".join(["_".join([str(j) for j in i]) for i in cov_std])
        # get avearge of intron covered length
        #infos = [[1,2,4], [3,6,9], [11,22,33]]
        #print "hhhhhhhhh", infos
        num_libs = len(infos[0])
        for j in range(num_libs):
            ave = numpy.mean([float(i[j]) for i in infos])
            #print ave
            all_lib_ave.append(ave)
        content =['all_libs'] + bb.split('_')
        writer.writerow(content)
        writer.writerow(["all_libs"]+ all_lib_ave)
        writer.writerow(["lib_name"]+exons_aa.split(':'))


def parse_args():
    parser = argparse.ArgumentParser(
        description='get coverage for a list of positions')
    parser.add_argument(
        '-m', '--meta_file',
        help='specify input file which is meta file',
        required=True)
    parser.add_argument(
        '-b', '--bed_regions',
        help='specify region of interest in bed format',
        required=True)
    # parser.add_argument(
    #     '-p', '--pairing', required=True,
    #     help='specify if sample paired with matched normal: '
    #           'paired or unpaired')
    args = parser.parse_args()
    return args



def __main__():
    # bam_files = "vcf_bam_all_samples_new_20150924.csv"
    # bed_regions = "WT1_combined_exon_pm_1kb.txt" 
    combined_file = "combined_exon_coverage.txt"
    exon_file = "WT1_intervals.txt"

    start = datetime.datetime.now()
    args = parse_args()
    meta_file = args.meta_file
    bed_regions =args.bed_regions
    input_headers = HEADER.INPUT_FILE_HEADER

    #first_exon_pos = 7565097
    first_exon_pos = 32409302
    patient_files = get_files(meta_file, input_headers)
    pprint(patient_files)


    print "Making BEDcoverage_scripts!\n"
    template_dir = '/home/szong/projects/development/coverage'
    
    out = make_BEDcoverage_script(patient_files, bed_regions, template_dir)
    BEDcoverage_scripts = out[0]
    complete_stamps = out[1]
    coverage_files = out[2]
    
    print "qsub BEDcoverage_scripts!\n"
    # qsub_scripts(BEDcoverage_scripts)

    print "Detect if cluster jobs finished!\n"
    # detect_cluster_jobs(complete_stamps) 
    header = ["chr", "start", "end", "exon", "pos"]
    #libs = [lib for lib in lib_bams]
    combine_coverage(coverage_files, combined_file, header)
    # sys.exit() 
    print "Make exon and intron dict from exon file!\n "
    out = make_exon_intron_dict(exon_file)
    exon_dict = out[0]
    intron_dict = out[1] 
    
    print "calculate per exon coverage!\n"
    print coverage_files,exon_dict, intron_dict, first_exon_pos
    calculate_avecov_all_libs(coverage_files, exon_dict, intron_dict, first_exon_pos)
    # sys.exit()

        
if __name__ == '__main__':
    __main__()

