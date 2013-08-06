#!/usr/bin/env python
from Bio import SeqIO
import random 
import re
import os
import sys

## Select random fasta sequence from fasta file. ##
seq_per_file = int(sys.argv[2])
record_dict  = SeqIO.index(sys.argv[1], "fasta")
selected_keys = random.sample(record_dict.keys(), seq_per_file)

for keys in selected_keys:
    handle = open("fasta_results_3.txt", "a")
    SeqIO.write(record_dict[keys], handle, "fasta")
    handle.close











