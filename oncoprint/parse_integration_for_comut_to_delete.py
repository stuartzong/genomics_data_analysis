#! /usr/bin/env python

import os, stat, os.path, time, datetime, subprocess
import re, sys, glob, argparse, csv, shutil, fileinput
from collections import defaultdict
from pprint import pprint
from itertools import islice
import ConfigParser

def __main__():
    print "script starts at: %s\n" % datetime.datetime.now()     
    parser = argparse.ArgumentParser(description='Parse integration results for comut plot!')
    parser.add_argument('-i','--input_file', help='specify integration summary file', required=True)
    args = vars(parser.parse_args())
    strain_list = "strains.txt"
    patient_list = "HIV_Cervical_patients.txt"
    integration_summary = args['input_file']
    comut_file = ".".join([integration_summary, "comut"])

    strains = make_list(strain_list)
    patients = make_list(patient_list)

    format_integrations(integration_summary, comut_file, strains, patients)


def make_list(infile):
    print "infile is: ", infile
    items = []
    with open(infile, 'r') as fh:
        for line in fh:
            item = line.split()[0]
            items.append(item)
        #items = list(set(items))
        print items
    return items

def format_integrations(integration_summary, comut_file, strains, patients):

    """ d = gene_integration_patients dictionary """
    d = dict()
    with open(comut_file,  'wb') as fh:
            headers = ['strains', 'sample', 'status']
            writer = csv.writer( fh, delimiter='\t' )
            writer.writerow(headers)
            print "The integration summary file is: %s.\n" % integration_summary
            with open(integration_summary, 'r') as handle:
                records = csv.DictReader(handle,  delimiter='\t')
                headers = records.fieldnames
                for line in records:
                    strain = "_".join(line['virus'].split())
                    patient = line["patient"].split()[0]
                    integration = line["Integration"]

                    try:
                        d[strain][patient].append(integration)
                    except KeyError:
                        if strain not in d:
                            d[strain] = {}
                        if patient not in d[strain]:
                            d[strain][patient] = [integration]
            print "yyyyyyyy"
            pprint(d)
            print "xxxxxx"
            for strain in strains:
                if (strain in strains):
                    int_status = []
                    num_patients = len(d[strain])
                    if (len(d[strain]) > 0):
                        for patient in patients:
                            #for patient in d[strain]:
                            if (patient in d[strain]):
                                status = list(set(d[strain][patient]))
                                #print status
                                int_status.append(status)
                                mut = list(set(d[strain][patient]))
                                if (int(len(mut)) > 1):
                                    integration = "Multiple"
                                else:
                                    integration = mut[0]
                            else:
                                integration = "NA"
                                mut = ["NA"]
                                #int_status = ["NA"]
                            #print strain, "\t",patient, "\t", integration,"\t", len(mut)
                            writer.writerow([strain, patient, integration])
                        final_status = [i[0] for i in int_status ]
                        integrated = final_status.count('YES')
                        unintegrated = final_status.count('NO')
                        print strain, "\t",num_patients, "\t", final_status, "\t", integrated, "\t", unintegrated
                        writer.writerow([strain, num_patients, integrated, unintegrated])


if __name__ == '__main__':
    __main__()

