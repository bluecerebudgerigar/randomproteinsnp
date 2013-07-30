#!/usr/bin/env python 


from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

import random 
import os
import subprocess
import re 
import getopt
import sys
import signal
from math import ceil




opts, arguments = getopt.getopt(sys.argv[1:], "f:", 
["fasta"])
for option, argument in opts:
    if option in ("-f", "--fasta"):
        fasta_file = argument
        record_dict    = SeqIO.index(argument, "fasta")
        
pwd = os.getcwd()

def find_cds ():
    seq_des = str(record_dict[keys].description).split("|")
    if any("CDS:" in s for s in seq_des):
        for des in seq_des:
            match = re.match("CDS:", des)
            if match is not None:
                print record_dict[keys].id
                feature, cds_start, cds_end = re.split(":|-", des)
                f = FeatureLocation(int(cds_start)-1, int(cds_end+3))
                cds_sequence = f.extract(record_dict[keys].seq)
                protein_sequence = cds_sequence.translate()
                if "*" not in protein_sequence:
                    return 3
                else:
                    return 1
        
    else:
        return 0
   

        
def write_file(object_name, file_name, mode):
    file_path=pwd + file_name
    handle = open(file_path, mode)    
    handle.write(object_name)
    handle.close()

i = 1
check = 0
for keys in record_dict:
    print i
    check = find_cds()
    fasta = ">%s\n%s\n" % (record_dict[keys].description, record_dict[keys].seq)
    
    if check == 0:
	write_file(fasta, "/error_nocds.txt", "a")
    elif check == 3:
        write_file(fasta, "/error_wrongcds.txt", "a")
    elif check == 1:
        write_file(fasta, "/cds_corrected_plus3.txt", "a")
    i +=1       
