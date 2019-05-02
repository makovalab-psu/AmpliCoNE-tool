#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd

parser=argparse.ArgumentParser(
    description='''Parse mapability scores and expand them to position specific values. ''')
parser.add_argument('--GC', '-g', required=True, help=' Chromosome wide GC-percentage values ', metavar='<TXT>')
parser.add_argument('--MAP', '-m', required=True, help=' Chromosome wide GEM-mappability values ', metavar='<TXT>')
parser.add_argument('--RM', '-r', required=True, help=' Chromosome wide REPEAT MASKED positions ', metavar='<TXT>')
parser.add_argument('--TRF', '-t', required=True, help=' Chromosome wide TANDEM REPEAT FINDER positions ', metavar='<TXT>')
parser.add_argument('--INFO', '-i', required=True, help=' Chromosome wide INFORMATIVE positions ', metavar='<TXT>')
parser.add_argument('--OUT', '-o', required=True, help=' Name of output file ', metavar='<NAME>')
parser.add_argument('--LENGTH','-l', type=int, required=True, help='Length of the Y chromosome in the reference.')
	
args=parser.parse_args()

GC=args.GC #Output from bin/calculate_Ychr_GCpercent_and_generate101bp_reads.py 
MAP=args.MAP #Output from bin/parse_Ychr_Mappability.py 
RM=args.RM #Output from parse_Ychr_RepeatMasker.py 
TRF=args.TRF #Output from bin/parse_Ychr_TRF.py 
INFO=args.INFO #Output from bin/parse_geneFamily_informativeSites.py 
len_chr=args.LENGTH #Length of Y chromosome in reference
OUT=args.OUT


# GC="/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/RM/chrY_GC_Content.txt"
# MAP="/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/GEM-MAPPABILITY/chrY_MAPPABILITY_bypos.txt"
# RM="/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/RM/chrY_REPEAT_MASKED.txt" 
# TRF="/nfs/brubeck.bx.psu.edu/scratch6/rahul/GTEx/REF/RM/chrY_TANDAMREPEAT_MASKED.txt"
# INFO="/nfs/brubeck.bx.psu.edu/scratch5/rahul/GTEx/TESTIS/REF/MAP_CUSTOM/Y_Informative_101bp_hg38.tab"
# len_chr=57227415

#read file as pandas dataframe

dfGC = pd.read_table(GC, sep="\t", header=0)
dfMAP = pd.read_table(MAP, sep="\t", header=0)
dfRM = pd.read_table(RM, sep="\t", header=0)
dfTRF = pd.read_table(TRF, sep="\t", header=0)
dfINFO = pd.read_table(INFO, sep="\t", header=0)

#range(1,len_chr+1) 1 - based output

df_Summary = pd.DataFrame()
df_Summary=df_Summary.assign(Position=range(1,len_chr+1))
df_Summary = pd.concat([df_Summary,dfGC,dfMAP,dfRM,dfTRF,dfINFO], axis=1, ignore_index=True)
df_Summary.columns=['Position','GCpercentage', 'Mappability', 'RepeatMask', 'TRF', 'InformativeSites']
df_Summary_norepeat=df_Summary[(df_Summary['RepeatMask']==1)]
df_Summary_norepeat=df_Summary_norepeat[(df_Summary_norepeat['TRF']==1)]
df_Summary_norepeat=df_Summary_norepeat[(df_Summary_norepeat['GCpercentage']>0.0)]
df_Summary_norepeat=df_Summary_norepeat[(df_Summary_norepeat['Mappability']>0.0)]
df_Annotation=df_Summary_norepeat.drop(['RepeatMask','TRF'],axis=1)
df_Annotation.to_csv(OUT, sep='\t', index=False)
