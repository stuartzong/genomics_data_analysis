#!/usr/bin/env python
import sys
import commands
import os
import argparse

__version__ = "1.0.2"

def create_job(job_file, cmds, output,name=None, email=None):
    #writes the job script
    out = open(job_file, 'w')
    out.write("#! /bin/bash\n")
    out.write("#$ -V\n")
    out.write("#$ -S /bin/bash\n")
    out.write("#$ -q thosts.q\n")
    out.write("#$ -j y\n")
    #out.write("#$ -P md5sum\n")
    out.write("#$ -o %s\n" % (output))

    if name:
        out.write("#$ -N %s\n" % (name))
    if email:
        out.write("#$ -M %s\n#$ -m e\n" % (email))

    out.write("set -euo pipefail\n")
    out.write('\n')
    out.write('echo "Job started at: $(date)"')
    out.write('\n')
    out.write('echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"')
    out.write('\n\n')

    for cmd in cmds:
        out.write(cmd + '\n')
    out.write('echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"')
    out.write('\n')
    out.write('echo "Job ended at:   $(date)"')
    out.write('\n')
    out.close()

    return True


def create_job_body(options, source, destination):
    cmds = []
    patient = source.split('/')[-1]
    cmds.append("echo " + '"' + "md5sum checking " + source + '"')
    #cmds.append("time rsync -vaHu --stats %s %s/" % (source, destination))
    cmds.append("find %s -type f | xargs md5sum > %s_md5sum_input" % (source, patient))
    cmds.append("echo " + '"' + "Successfully generated md5sum for" + source + "\n")
    return cmds


def check_script_created(successful):
    if successful:
        print "Wrote the script: %s.sh" % (job_name)
    else:
        print "Error in writing job script!"


#Options for script
parser = argparse.ArgumentParser(epilog="Tips: Run script in the directory you want the log files to be. Be sure to include absolute paths to -n and -o.",
                                 version= "transfer " + __version__)
required_mutually_exclusive = parser.add_argument_group('required mutually exclusive arguments')
mutually_exclusive = required_mutually_exclusive.add_mutually_exclusive_group(required=True)
mutually_exclusive.add_argument("-n", dest="input",
                                help="directory or file to transfer; only accepts a single file or directory", metavar="SOURCE")
mutually_exclusive.add_argument("-N", dest="inputs", nargs="+",
                                help="files to transfer; accepts multiple files", metavar="SOURCES")
parser.add_argument("-o", dest="output_directory",
                   help="destination directory", metavar="DESTINATION")
parser.add_argument("-e", "--email", dest="mail",
                   help="address to send email upon completion", metavar="EMAIL") #mail is the email address
parser.add_argument("-l", dest="directory",
                   action="store_true", default=False,
                   help="transfer the contents of the given input directory separately instead of the entire directory at once (creates transfer script for each item in the given input directory)")
parser.add_argument("-p", dest="permissions",
                   action="store_true", default=False,
                   help="check permissions to see if files could be removed")
parser.add_argument("-i", dest="id", default="",
                   help="id to help make the transfer script file name unique")
options = parser.parse_args()

id = options.id
admin = None

if id != "":
    id = id + "-"
transfertype = "md5sum"

#needs input and output directory
if options.input and os.path.isfile(options.input) and options.directory:
    print "Error, can't use the -l option if a file is passed to the -n option since a file can't be split into multiple transfers."
    sys.exit()
elif options.input and options.output_directory:
    source = os.path.abspath(options.input)
    basename = os.path.basename(source)
    dest = os.path.abspath(options.output_directory)
    destination = os.path.join(dest,basename)
    if not os.path.exists(source):
        print "input path does not exist"
        sys.exit()
    if os.path.exists(destination):
        print 'Warning, a directory "' + basename + '" already exists in ' + dest
    if not os.path.exists(dest):
        print 'Warning, ' + dest + ' does not already exist, it will be created with "' + basename + '" inside it'
    if options.mail:
        admin = options.mail
elif options.inputs and options.output_directory:
    for source in options.inputs:
        if os.path.isfile(source):
            source_abs_path = os.path.abspath(source)
            basename = os.path.basename(source_abs_path)
            dest = os.path.abspath(options.output_directory)
            destination = os.path.join(dest,basename)
            if not os.path.exists(source):
                print "input path does not exist"
                sys.exit()
            if os.path.exists(destination):
                print 'Warning, a directory "' + basename + '" already exists in ' + dest
            if not os.path.exists(dest):
                print 'Warning, ' + dest + ' does not already exist, it will be created with "' + basename + '" inside it'
            if options.mail:
                admin = options.mail

            #transfer the entire given input directory
            cmds = create_job_body(options, source_abs_path, dest)
            job_name = os.path.basename(source_abs_path) + "-" + id + transfertype
            return_status = create_job(job_name + ".sh", cmds, os.getcwd(), name=job_name, email=admin)
            check_script_created(return_status)
        else:
            print("ERROR: -N only works with files, given: " + source)
    sys.exit()
else:
    print "Need input(-n), destination(-o). Use the -h option for help"
    sys.exit()

if options.directory:
    # split the contents of the given input directory into separate transfer scripts
    for dirs in os.listdir(source):
        pathname = os.path.join(source, dirs)
        parent_directory = os.path.basename(source)
        modified_dest = os.path.join(dest, parent_directory)
        cmds = create_job_body(options, pathname, modified_dest)
        job_name = parent_directory + "-" + dirs + "-" + id + transfertype
        return_status = create_job(job_name + ".sh", cmds, os.getcwd(), name=job_name, email=admin)
        check_script_created(return_status)
else:
    #transfer the entire given input directory
    cmds = create_job_body(options, source, dest)
    job_name = os.path.basename(source) + "-" + id + transfertype
    return_status = create_job(job_name + ".sh", cmds, os.getcwd(), name=job_name, email=admin)
    check_script_created(return_status)

if options.permissions:
    print "Checking permissions of original files..."
    output = commands.getstatusoutput("python /home/cchoo/scripts/check_permissions5.py -n %s" % (source))
    print output[1]
