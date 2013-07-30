#!/usr/bin/env python

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import random 
import os
import subprocess
import re 
import getopt
import sys
import signal
from math import ceil

from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

opts, arguments = getopt.getopt(sys.argv[1:], "f:i:c:", 
["fasta","iterations","cutoff"])
for option, argument in opts:
    if option in ("-f", "--fasta"):
        fasta_file = argument
        record_dict    = SeqIO.index(argument, "fasta")
    elif option in ("-i", "--iterations"):
        iterations     = int(argument)
    elif option in ("-c", "--cutoff"):
        transcripts_cutoff = float(argument)

        
pwd = os.getcwd()
print fasta_file

def check_syn (ori_AA, mut_AA, base_loc):
    n = 0 
    aa_input = ""
    if ori_AA == '*' or mut_AA == '*':
        n = 2
    elif ori_AA != mut_AA:
        base_loc = float(base_loc)
        aa_pos = base_loc/3
        aa_pos = ceil(aa_pos)
        aa_pos = int(aa_pos)
        print aa_pos
        aa_input = "%s%s%s" % (ori_AA, str(aa_pos), mut_AA) 
        n = 1
    elif ori_AA ==  mut_AA:
        n = 0
    
    return n, aa_input
        
def snp_translation (snp_base, base_loc):
    i = int(base_loc) ##i = 124
    if (i%3 == 0):
        original_seq = str(cds_sequence)[i-3:i]
        mutant_seq = original_seq[:2] + snp_base
    elif (i%3 == 1): ## FUCKING CHANGE HERE NEED TO SHIT FIRST 
        original_seq = str(cds_sequence)[i-1:i+2] ##126
        mutant_seq = snp_base + original_seq[1:3]
    elif (i%3 == 2):
        original_seq = str(cds_sequence)[i-2:i+1]
        mutant_seq = original_seq[0] + snp_base + original_seq[2]
        
    print "ori = %s mut = %s" % (original_seq, mutant_seq)
    ori_AA = str(Seq(original_seq, IUPAC.unambiguous_dna).translate())
    mut_AA = str(Seq(mutant_seq, IUPAC.unambiguous_dna).translate())
    return ori_AA, mut_AA
    
def find_cds ():
    seq_des = str(record_dict[keys].description).split("|")
    for i in seq_des:
        if re.match("CDS", i):
            feature, cds_start, cds_end = re.split(":|-", i)
    cds_feature = SeqFeature(FeatureLocation(int(cds_start)-1,int(cds_end)+1),
                type=str(feature))
    cds_sequence = cds_feature.extract(record_dict[keys].seq)
    print cds_sequence
    protein_sequence = cds_sequence.translate()
    return cds_start, cds_end, cds_sequence, protein_sequence

def prepare_command_line ():
    seq_name = record_dict[keys].id
    supporting_set = seq_name + ".sss"
    provean_cmd = ["provean.sh","-q","fasta_record","-v","var_file_2"]
    if os.path.isfile(supporting_set):
        provean_cmd.append("--supporting_set")
    else:
        provean_cmd.append("--save_supporting_set")
    provean_cmd.append(supporting_set)
    return provean_cmd
    
def random_snps(cycle_number):
    snp_input_list = []
    snp_input      = ""
    snp_pos=random.sample(range(1,cds_length+1), cycle_number)
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
neg_hits = 0 



for keys in record_dict:
    non_hits = 0
    cds_start, cds_end, cds_sequence, protein_sequence = find_cds()
    print "first base = " + cds_sequence[0]
    print "first aa = " + protein_sequence[0]
    cds_length = len(str(cds_sequence))
    max_snps = cds_length * 3 
    cds_name = record_dict[keys].id
    if max_snps >= iterations:
        neg_hits_cutoff = iterations - transcripts_cutoff
    elif max_snps < iterations:
        transcripts_cutoff = int(ceil(0.95 * float(max_snps)))
        neg_hits_cutoff = max_snps - transcripts_cutoff
        
    print cds_sequence
    fasta_record = ">%s\n%s" % (cds_name, protein_sequence)
    print fasta_record
    write_file(fasta_record, "/fasta_record", 'w')
    i2 = 0
    i = 1
    
    
    while i2 < 1:  ###(master control for sequences)
        n        = 0
        hits     = 0
        true_neg_hits = 0
        pos_hits = 0
        
        print "checkpoint0"
        snp_list = []
        while (n < iterations) and (n < cds_length*3) : ###(iternations and number of snps to introduced.)
            print "checkpoint1"
            print len(snp_list)
            print cds_length 
            snp_input = random_snps(i)
            print snp_input
            if snp_input not in snp_list:
                snp_list.append(snp_input)
                n = n + 1
            else:
                continue
        n = 0 # Reset n counter for the next flow control  
        for snps in snp_list: #take each snp, and divde them, ie. if 2snp/run = 
            print "checkpoint2"
            print snps
            print cds_length
            snps_input = snps.split("-")
            no_of_snps = len(snps_input)
            master_control = 0
            for snps in snps_input: #2nd level of snp, if any of the snp is del then pos_hits + 1
                print "checkpoint #5"
                snp_base = snps[-1]   # get the mutant base
                base_loc = snps[1:-1] # get the base location
                print "base loc = %s snp_base == %s ori_base = %s" % ( str(base_loc), snp_base, snps[0])
                ori_AA, mut_AA = snp_translation(snp_base, base_loc) 
                print "ori = %s, mut = %s" % (ori_AA, mut_AA)
                n, aa_input = check_syn(ori_AA, mut_AA, base_loc)
                print aa_input
                print n
                if (base_loc < 4) or (n == 2):
                    master_control = 1
                    break
                elif n == 1:
                    write_file(aa_input + "\n", "/var_file_2", "w")
                    provean_cmd=prepare_command_line()
                    results = run_provean()
                    print results
                    if float(results[-1]) < (-2.5) :
                        master_control = 1
                        break
            if master_control != 1:
                true_neg_hits += 1   
                if true_neg_hits > neg_hits_cutoff:
                    break
            if master_control == 1:
                pos_hits = pos_hits + 1
                if pos_hits == transcripts_cutoff: 
                    end_results = "%s\t%s\t%s\n" % (cds_name, str(i), cds_length)
                    write_file(end_results, "/results", "a")
                    i2 = 1
                    break
            print "pos hits = %s, neg_hits = %s" % (str(pos_hits), str(true_neg_hits))
        i += 1
        
            

            


            
        








