#!/usr/bin/python

# written by Vipin T Sreedharan, Gunnar Raetsch and Andre Kahles

# Usage : 
# 1. Provide a valid bed file.
# 2. Number of lines to be extracted from the file.
# 3. Provide result file in bed format.

import sys
import random 

# Accept a bed file from user
in_file = sys.argv[1]
# Number of lines to be extracted from the bed file 
target_no = int(sys.argv[2])
# Result file 
out_file = sys.argv[3]

print("")
print("subsample_bed.py version 0.1")
print("------------------------------")

# open the file 
no_lines = 0 
for line in open(in_file, "rU"):
    no_lines += 1

try:
    accept_prob = (1.0 * target_no) / no_lines
except:
    print("\n  Computing acceptance probability failed, setting to 1\n") 
    accept_prob = 1 

# Input verification 
if target_no > no_lines:
    print("\n  Warning:\tLess lines present in bed file than target number of lines to extract.")
    print("\t\tTotal number of lines in fastq file : " + str(no_lines))
    print("\t\tNumber of lines to be extracted : " + str(target_no) + "\n")
    accept_prob = 1 
	
# Taking random lines.    
print "  Total number of lines present in the input file : ", no_lines
print "  Probability of accepting a line:                        ", accept_prob
res_line = 0
c = 0 
res_file = open(out_file, "w+")
#for rec in SeqIO.parse(in_handle, "fastq"):
for line in open(in_file, "rU"):
    c += 1 
    if random.random() <=  accept_prob:
        res_line += 1
        print >> res_file, line,

    if target_no == res_line: # stop when we have enough lines
        break
print "  Number of lines taken :                                 ", res_line
print "  Total number of lines scanned :                         ", c
#print '----------------------------------------------'
res_file.close()
