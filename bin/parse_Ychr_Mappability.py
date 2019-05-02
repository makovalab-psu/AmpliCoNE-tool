#!/usr/bin/python
import sys
import argparse
import os.path
import pandas as pd
parser=argparse.ArgumentParser(
    description='''Parse mapability scores and expand them to position specific values. ''')
parser.add_argument('--INPUT', '-i', required=True, help=' GEM-mappability output in bed format ', metavar='<BED>')
parser.add_argument('--CHR', '-c', required=True, help=' Chromosome annotation as found in the GEM-mappability output file. ', choices=['Y','chrY'])
parser.add_argument('--LENGTH','-l', type=int, required=True, help='Length of the Y chromosome in the reference.')
	
args=parser.parse_args()

INPUT=args.INPUT #GEM-mappability output in bed format
CHR=args.CHR #chromosome name as defined in the mappability output Y or chrY
len_chr=args.LENGTH #Length of Y chromosome in reference

#read file as pandas dataframe
df = pd.read_table(INPUT, sep="\t", header=None, low_memory=False)
chrY=df[(df[0]==CHR)] #parse Y chromsome specific reads
del [df] #Free up space

#Update column to create proper intervals in each row which represents the window with mappability value.
chrY[1]+=1 
temp=chrY[2].iloc[1:].tolist()
temp.append(len_chr+1)
chrY[2]=temp

#Expand the windows to individual positions
map_val = []
for row in chrY.itertuples(index=False): #By setting the index parameter to False we can remove the index as the first element of the tuple.
	map_val+=[float(row[4])]*(int(row[2])-int(row[1]))

#print the file
file = open(str(CHR)+"_MAPPABILITY_bypos.txt", "w")
file.write("Mappability\n")
for index in range(len(map_val)):
    file.write(str(map_val[index])+"\n")


file.close()
