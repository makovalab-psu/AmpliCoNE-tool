#!/usr/bin/python
import sys
import argparse

parser=argparse.ArgumentParser(
    description='''Parse REPEAT MASKER annotation and expand them to position specific values on Y. ''')
parser.add_argument('--INPUT', '-i', required=True, help=' REPEAT MASKER output format ', metavar='<RM-OUT>')
parser.add_argument('--CHR', '-c', required=True, help=' Chromosome annotation as found in the REPEAT MASKER output file. ', choices=['Y','chrY'])
parser.add_argument('--LENGTH','-l', type=int, required=True, help='Length of the Y chromosome in the reference.')
	
args=parser.parse_args()

INPUT=args.INPUT #REPEAT MASKER output
CHR=args.CHR #chromosome name as defined in the mappability output Y or chrY
len_chr=args.LENGTH #Length of Y chromosome in reference

#INPUT='/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/RM/hg38.sorted.fa.out'

repeatmask_val = [1]*len_chr
with open(INPUT,"r") as rFH:
	next(rFH)#heading
	next(rFH)#heading
	next(rFH)#empty line
	for line in rFH :
		col=line.rstrip('\n').split()
		if col[4] == CHR:
			repeatmask_val[int(col[5]):int(col[6])]=[0]*(int(col[6])-int(col[5]))

rFH.close()

file = open(str(CHR)+"_REPEAT_MASKED.txt", "w")
file.write("RepeatMasked\n")
for index in range(len(repeatmask_val)):
	file.write(str(repeatmask_val[index])+ "\n")

file.close()
