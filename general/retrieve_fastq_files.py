#!/projects/trans_scratch/software/python/Python-2.7.3/bin/python

# Script to retrieve the merged BAM file for the given library name.

import argparse, glob, os, subprocess, sys, xmlrpclib
from gsc import getr2gbam

try:
    USERNAME = os.environ['USER']
except KeyError:
    print("ERROR: Cannot find USER in the environment!")
    exit(1)

try:
    PACKAGE_DIRECTORY = os.environ['TRANSABYSS_PATH']
except KeyError:
    print("ERROR: Cannot find TRANSABYSS_PATH in the environment!")
    exit(1)

BIOAPPS_API_SERVER = 'www.bcgsc.ca/data/sbs/viewer/api'

FAILED_MESSAGE = 'Failed to find anything =('
REFERENCE_CHECKER = '/projects/trans_scratch/validations/assembly_scripts/misc/get_reference_from_bam.sh'

def searchForBAMFiles(directory):
    """
    Search for BAM files in the given directory.

    @param directory - the directory to search for BAM files.

    @return True if at least one BAM file is found in the given directory.

    """
    found = False
    bam_files = glob.glob(os.path.join(directory, "*bam"))
    for bam_file in bam_files:
        bashCommand = REFERENCE_CHECKER + " " + bam_file
        process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
        output = process.communicate()[0].rstrip()
        print(output + " BAM file found: " + bam_file)
        index_file = bam_file + ".bai"
        if not "incomplete.bam" in index_file:
            if os.path.isfile(index_file):
                print("BAM index found: " + index_file)
            else:
                print("Missing BAM index file for: " + bam_file)
        found = True
    return found

def retrieveFastqFiles(api_fastq_result, production_only=False):
    target_number_of_lanes=api_fastq_result[0]['target_number_of_lanes']
    if target_number_of_lanes == None:
        return None
    data_paths = []
    for lane_info in api_fastq_result:
        for fastq in lane_info['fastq']:
            if production_only:
                if fastq['status'] != 'production':
                    continue
            data_path=fastq['data_path']
            if "chastity_passed" in data_path:
                data_paths.append(data_path)
    fastq_files=[]
    for path in data_paths:
        path_found = searchForFASTQFiles(path)
        if path_found == True:
            fastq_files.append(path)
    expected = target_number_of_lanes*2
    if len(fastq_files) > expected:
        print("WARNING: MORE FASTQ FILES WERE FOUND ON FILE SYSTEM THAN EXPECTED")
    elif len(fastq_files) < expected:
        print("WARNING: LESS FASTQ FILES WERE FOUND ON FILE SYSTEM THAN EXPECTED")
    else:
        print("ALL FASTQ FILES FOUND:")
    for i in fastq_files:
        print(i)
    return fastq_files

def searchForFASTQFiles(path):
    """
    Search for an existing fastq file given a file path

    @param path - the path to check for a fastq file.

    @return True if the FASTQ file is found.
    """
    return True if os.path.isfile(path) and os.path.getsize(path) > 0 else False


def get_fastq_read_lengths(api_fastq_result):
    sequence_len = []
    fastq_info = [i['fastq'] for i in api_fastq_result]
    for i in fastq_info:
        if len(i) == 0:
            continue
        else:
            for j in i:
                sequence_len.append(j['sequence_length'])
    return sequence_len

parser = argparse.ArgumentParser(description='Retrieve the filer path of the merged or repositioned BAM file for the given library name.',
                                 add_help=False)
required = parser.add_argument_group('POSITIONAL ARGUMENTS')
required.add_argument('library_name', metavar='LIBRARY_NAME', help='Library name of the merged or repositioned BAM file to retrieve, e.g. A10538.')
optional = parser.add_argument_group('OPTIONAL ARGUMENTS')
optional.add_argument('-h', '--help', action='help', help='show this help message and exit')
optional.add_argument('-q', '--qc',action='store_true', default='False', help='Do not Filter the fastq files that have not passed qc.')
optional.add_argument('-f', '--filter', action='store_true', default=False, help='Do not Filter the fastq files that are not marked as production.')
if len(sys.argv) == 1:
    parser.print_help()
    exit(2)

args = parser.parse_args()
library_name = args.library_name
ignore_qc = args.qc
production_only = not args.filter
try:
    bioapps_api = xmlrpclib.ServerProxy("http://%s:%s@%s" % (USERNAME, getr2gbam.getpassword(), BIOAPPS_API_SERVER), allow_none=True)
    api_fastq_results = bioapps_api.getFASTQInfo({'library':library_name,'ignore_qc':ignore_qc})
    if api_fastq_results:
        read_lengths = get_fastq_read_lengths(api_fastq_results)
        if len(read_lengths) == 0:
            print(FAILED_MESSAGE)
        for read_length in set(read_lengths):
            print("Retreiving fastqs for read length {}".format(read_length))
            filtered_api_fastq_result = bioapps_api.getFASTQInfo({'library':library_name,'ignore_qc':ignore_qc,'sequence_length':read_length})
            fastqs = retrieveFastqFiles(filtered_api_fastq_result,production_only)
    else:
        print(FAILED_MESSAGE)
except xmlrpclib.ProtocolError as error:
    print('ERROR: There appears to be an issue accessing the BIOAPPS database, "{0} {1}"'.format(error.errcode, error.errmsg))
    exit(1)
