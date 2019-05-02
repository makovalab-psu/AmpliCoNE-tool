#!/usr/bin/python
import sys
import argparse
from Bio import SeqIO


parser=argparse.ArgumentParser(
    description='''Caliculate GC percentage of Y chromosome using a sliding window of 250 ''')
parser.add_argument('--REF', '-r', required=True, help=' Reference genome', metavar='<FASTA>')
parser.add_argument('--CHR', '-c', required=True, help=' Name of chromosome as found in FASTA header. ', choices=['Y','chrY'])

	
args=parser.parse_args()

INPUT=args.REF #Reference genome
CHR=args.CHR #chromosome name as defined in the Reference genome Y or chrY


def chunks(seq, window):
	seqlen = len(seq)
	for i in range(0,seqlen,window):
		if (i+window>seqlen):
			j = seqlen 
		else:
			j = i+window
		if i-250 < 0:
			l=0
		else:
			l=i-250
		if j+250 > seqlen:
			r=seqlen
		else:
			r=j+250
		yield seq[l:r],i+1,j
		if j==seqlen:
			break


#INPUT='/galaxy/home/rxv923/refs/Hg/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
refFH=open(INPUT,"rU")
chr_len=0
for chr_record in SeqIO.parse(refFH, "fasta"):
	if chr_record.id==CHR:
		chr = chr_record.seq
		chr_len=len(chr)
		GC=[]
		for win in chunks(chr, 1):
			N = win[0].count("N")
			A = win[0].count("A")
			T = win[0].count("T")
			G = win[0].count("G")
			C = win[0].count("C")
			#id=str(win[1])+"-"+str(win[2])
			GC+=[((G+C)/float(A+T+G+C+0.0000001))*100]
		

refFH.close()

file = open(str(CHR)+"_GC_Content.txt", "w")
file.write("GC\n")
for index in range(len(GC)):
	file.write(str(GC[index])+ "\n")

file.close()

file2 = open(str(CHR)+"_length.txt", "w")
file2.write(str(chr_len)+ "\n")
file2.close()


#read length to generate reads for informative sites calculation
r_len=101

fasta='Sliding_window101bp_Reference_reads.fasta'
#Generate fastq file with each read header has chromosome name and position information (>chr_position) and a sequence is 101bp which starts at position defined in header.
fasta_output=open(fasta,"w")
for i in range(len(chr)-r_len):
	r_seq=str(chr[i:i+r_len])
	r_header=">"+str(CHR)+"_"+str(i+1)
	fasta_output.write(r_header+"\n")
	fasta_output.write(r_seq+"\n")

fasta_output.close()
