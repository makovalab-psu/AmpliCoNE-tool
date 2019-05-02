#!/usr/bin/python
import sys
import argparse

parser=argparse.ArgumentParser(
    description='''Parse TANDEM REPEAT FINDER annotation and expand them to position specific values on Y. ''')
parser.add_argument('--INPUT', '-i', required=True, help=' TANDEM REPEAT FINDER output', metavar='<BED>')
parser.add_argument('--CHR', '-c', required=True, help=' Chromosome annotation as found in the REPEAT MASKER output file. ', choices=['Y','chrY'])
parser.add_argument('--LENGTH','-l', type=int, required=True, help='Length of the Y chromosome in the reference.')
	
args=parser.parse_args()

INPUT=args.INPUT #TANDEM REPEAT FINDE output
CHR=args.CHR #chromosome name as defined in the mappability output Y or chrY
len_chr=args.LENGTH #Length of Y chromosome in reference

#INPUT='/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/RM/trfMask.bed'

trf_val = [1]*len_chr
with open(INPUT,"r") as trFH:
	for line in trFH :
		col=line.rstrip('\n').split('\t')
		if col[0] == CHR:
			trf_val[int(col[1]):int(col[2])]=[0]*(int(col[2])-int(col[1]))


trFH.close()

file = open(str(CHR)+"_TANDAMREPEAT_MASKED.txt", "w")
file.write("TRF\n")
for index in range(len(trf_val)):
	file.write(str(trf_val[index])+ "\n")

file.close()
