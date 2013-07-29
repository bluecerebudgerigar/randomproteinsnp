#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import random 
import os
import subprocess
import re 
import getopt

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

opts, arguments = getopt.getopt(sys.argv[1:], "f:i:c:", 
["fasta","iterations","pvalue","cutoff","quick"])
for option, argument in opts:
    if option in ("-f", "--fasta"):
        record_dict    = SeqIO.index(argument, "fasta")
    elif option in ("-i", "--iterations"):
        iterations     = int(argument)
    elif option in ("-c", "--cutoff"):
        transcripts_cutoff = float(argument)

neg_hits_cutoff = iterations - transcripts_cutoff
        


def check_syn (ori_AA, mut_AA):
    if ori_AA == '*' or mut_AA == '*':
      n = 2
    elif ori_AA != mut_AA:
        aa_input = "%s%s%s" % (ori_AA, str(i), mut_AA) 
        n =+ 1
    
    return n, aa_input
        
def snp_translation (snp_base, base_loc):
    i = int(base_loc)
    if (i%3 == 0):
        original_seq = str(cds_sequence)[-3:i]
        mutant_seq = original_seq[:2] + snp_base
    elif (i%3 == 1):
        original_seq = str(cds_sequence)[i-1:i+2]
        mutant_seq = original_seq[0] + snp_base + original_seq[2]
    elif (i%3 == 2):
        original_seq = str(cds_sequence)[i-2:i+1]
        mutant_seq = snp_base + original_seq[1:3]
        
    ori_AA = str(Seq(query_sequence, IUPAC.unambiguous_dna).translate())
    mut_AA = str(Seq(mutant_seq, IUPAC.unambiguous_dna).translate())
    return ori_AA, mut_AA
    
def find_cds ():
    seq_des = str(record_dict[keys].description).split("|")
    for i in seq_des:
        if re.match("CDS", i):
            feature, cds_start, cds_end = re.split(":|-", i)
    cds_feature = SeqFeature(FeatureLocation(int(cds_start)-1,int(cds_end)-1),
                type=str(feature))
    cds_sequence = cds_feature.extract(record_dict[keys].seq)
    print cds_sequence.translate()
    return cds_start, cds_end, cds_sequence

def prepare_command_line ():
    seq_name = record_dict[keys].id
    supporting_set = seq_name + ".sss"
    provean_cmd = ["provean.sh","-f","fasta_1.txt","-v",var_file]
    if os.path.isfile(supporting_set):
        provean_cmd.append("--supporting_set")
    else:
        provean_cmd.append("--save_supporting_set")
    provean_cmd.append(supporting_set)
    return provean_cmd
    
def random_snps(cycle_number):
    snp_input_list = []
    snp_input      = ""
    snp_pos=random.sample(range(1,cds_length), cycle_number)
    snp_pos.sort()
    for i in snp_pos:
        base      = cds_sequence[i-1]
        snp_base  = random.choice([x for x in "ATGC" if x != base])
        snp_input = "%s%s%s" % (base, str(i), snp_base)
        snp_input_list.append(snp_input)
    snp_input="-".join(snp_input_list)
    return snp_input

def run_provean():
    proc = subprocess.Popen(provean_cmd, stdout=subprocess.PIPE)
    results = proc.stdout.read()
    results = results.split("\n")
    results = results[-2].split("\t")
    return results
    
def write_file(object_name, file_name, mode):
    file_path=pwd + file_name
    handle = open(file_path, mode)    
    handle.write(object_name)
    handle.close()
    
base = ""
snp_base = ""
cutoff = 5
cds_start = 0


for keys in record_dict:
    non_hits = 0
    cds_start, cds_end, cds_sequence = find_cds()
    cds_length = len(str(cds_sequence))
    cds_name = record_dict[keys].id
    
    
    
    while i2 < 1:  ###(master control for sequences)
        n        = 0
        hits     = 0
        pos_hits = 0
        snp_list = []
        while n < iterations: ###(iternations and number of snps to introduced.)
            snp_input = random_snps(i)
            if snp_input not in snp_list:
                snp_list.append(snp_input) 
                n = n + 1
            else:
                continue
        n = 0 # Reset n counter for the next flow control  
        for snps in snp_list: #take each snp, and divde them, ie. if 2snp/run = 
            some_other_control = 0    
            snps_input = snps.split("-")
            no_of_snps = len(snps_input)
            master_control = 0
            while master_control < no_of_snps :
                for snps in snp_input: #2nd level of snp, if any of the snp is del then pos_hits + 1
                    snps.pop(0) # removes the orignal base.
                    snp_base = snps.pop(-1) # get the mutant base
                    base_loc = snp # get the base location
                    if base_loc < 4:
                        some_other_control = 1
                        pos_hits =+ 1
                        master_control = no_of_snps
                        break
                    else :
                        ori_AA, mut_AA = snp(translation, snp_base, base_loc) 
                        n = check_syn(ori_AA, mut_AA)
                        if n == 2:
                            pos_hits =+1
                            master_control = no_of_snps
                            break
                        elif n == 1:
                            write_file(aa_input, var_file, "w")
                            some_other_contorl = 1
                            provean_cmd=prepare_command_line()
                            results = run_provean()
                            if int(results[-1]) < -2.5 :
                                pos_hits =+ 1
                                master_control = no_of_snp
                if n == 0:
                    neg_hits = +1     
                if neg_hits > neg_hits_cutoff: 
                    i =+ 1
                if pos_hits == transcripts_cutoff: 
                    i2 = 1 
        write_file(i, results, "a")

            


            
        








