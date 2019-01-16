import os
import re
import argparse
import subprocess

from os.path import join

def walk_files(wkdir, old, new):
#     walk the dir tree down to top
    for root, dirs, files in os.walk(wkdir, topdown=False):
        for fname in files:
            fpath = join(root, fname)
            change_file_content(fpath, old, new)
            searchobj = re.search(old, fname, flags=0)
            if searchobj:
                newfname = fname.replace(old, new)
                newfpath = join(root, newfname)
                print "rename file: %s to \n%s" % (fpath, newfpath)
                os.rename(fpath, newfpath)
        for dname in dirs:
            dpath = join(root, dname)
            searchobj = re.search(old, dname, flags=0)
            if searchobj:
                newdname = dname.replace(old, new)
                newdpath = join(root, newdname)
                print "rename directory: %s to \n%s" % (dpath, newdpath)
                os.rename(dpath, newdpath)
            print dpath

def change_file_content(afile, old, new):
    # print "change file content using sed command for\n%s" % afile
    p = subprocess.Popen('sed -i -s "s/%s/%s/g" %s' % (old, new, afile),
                         shell=True, stdout=subprocess.PIPE)
    output,  err = p.communicate()


def parse_args():
    parser = argparse.ArgumentParser(
        description='change file conent, file name and dir name based on a map! ')
    parser.add_argument(
        '-i', '--input_file',
        help='specify input file with 3 columns: dir, old name, new name!',
        required=True)
    args = parser.parse_args()
    return args


def __main__():
    args = parse_args()
    infile = args.input_file
    # wkdir = "/projects/trans_scratch/validations/genome-validator/NCI_ALL_MPAL/test/TARGET-10-SJMPAL043771_gv2/gv-2.3"
    with open(infile, 'r') as handle:
        for line in handle:
            wkdir, old, new = line.strip().split('\t')
            print "start renaming dir: %s" % wkdir
            walk_files(wkdir, old, new)


if __name__ == '__main__':
    __main__()
