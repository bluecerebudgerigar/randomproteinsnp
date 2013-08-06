#!/usr/bin/env python
from Bio import SeqIO
import random 
import re
import os
import sys

## Select random fasta sequence from fasta file. ##
seq_per_file = int(sys.argv[2])
record_dict  = SeqIO.index(sys.argv[1], "fasta")
number_of_seq =int(sys.argv[3])
selected_keys = random.sample(record_dict.keys(), number_of_seq)


n = 0 
file_number = 1
for keys in selected_keys:
    if n < seq_per_file:
        fasta_file_name = "fasta_results_%i.fa" % (file_number)
        handle = open(fasta_file_name, "a")
        SeqIO.write(record_dict[keys], handle, "fasta")
        handle.close
        n += 1
    elif n >= seq_per_file:
        n = 0
        file_number += 1
        fasta_file_name = "fasta_results_%i.fa" % (file_number)
        handle = open(fasta_file_name, "a")
        SeqIO.write(record_dict[keys], handle, "fasta")
        handle.close
        n += 1
        
    
        
        
    










