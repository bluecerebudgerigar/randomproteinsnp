#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import random 

def translation (query_sequence, mutant_seq):
    global n
    global aa_input
    global del_snp 
    ori_AA = str(Seq(query_sequence, IUPAC.unambiguous_dna).translate())
    mut_AA = str(Seq(mutant_seq, IUPAC.unambiguous_dna).translate())
    if ori_AA == '*' or mut_AA == '*':
        del_snp =+ 1
    elif ori_AA != mut_AA:
        aa_input = "%s%s%s" % (ori_AA, str(i), mut_AA) 
        n =+ 2
def check_syn (snp_base):
    if (i%3 == 0):
        original_seq = str(sequence)[i-3:i]
        mutant_seq = original_seq[:2] + snp_base
        translation(original_seq, mutant_seq)
    elif (i%3 == 1):
        original_seq = str(sequence)[i-1:i+2]
        mutant_seq = original_seq[0] + snp_base + original_seq[2]
        translation(original_seq, mutant_seq)
    elif (i%3 == 2):
        original_seq = str(sequence)[i-2:i+1]
        mutant_seq = snp_base + original_seq[1:3]
        translation(original_seq, mutant_seq)
    return n
    return aa_input
    
base = ""
snp_base = ""
record_dict    = SeqIO.index("fasta_1.txt", "fasta")    

for keys in record_dict:
    sequence = record_dict[keys].seq 
    snp_pos=random.sample(range(1,len(str(sequence))), 10)
    snp_pos.sort()
    print snp_pos
    
    aa_input_list = []
    for i in snp_pos:
        aa_input=""
        n = 0
        del_snp = 0
        base = str(sequence)[i-1]
        snp_base = random.choice([x for x in "ATGC" if x != base])
        if check_syn(snp_base) == 2:
            aa_input_list.append(aa_input)
        
    print aa_input_list
    print del_snp
    
    


            


            
        








