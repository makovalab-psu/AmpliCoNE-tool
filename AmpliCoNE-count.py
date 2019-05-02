#!/usr/bin/python
import sys
import re
import pandas
import numpy
import subprocess
import copy
import argparse
import os
import pysam

parser=argparse.ArgumentParser(
    description='''Ampliconic Gene Copy Number Estimator (AmpliCoNE-copy) : Estimates the copy number of the 9 ampliconic gene families on Human Y chromosome. ''',
    epilog="""Email: v.rahul.simham@gmail.com for errors and bugs""")
parser.add_argument('--BAM', '-b', required=True, help=' Indexed BAM file. (Aligned using BWA MEM) ', metavar='<BAM>')
parser.add_argument('--CHR', '-c', required=True, help=' Chromosome annotation as found in the BAM header file. ', choices=['Y','chrY'])
parser.add_argument('--GENE_DEF', '-g', required=True, help=' Gene family and control gene definition ')
parser.add_argument('--ANNOTATION', '-a', required=True, help=' Y Chromosome annotation (GCper, Mappability,InformativeSites) ')
parser.add_argument('--LENGTH','-l', nargs='?', type=int, default=57227415, help='Length of the Y chromosome in the reference (hg38).(default: %(default)s) ')
parser.add_argument('--OUTPUT','-o', nargs='?', default='Output', help='Length of the Y chromosome in the reference (hg38).(default: %(default)s) ')
parser.add_argument('--READ','-r', nargs='?', default="PAIRED", help='The reads are paired end or single end, if paired we filter for proper read pairs. (default: %(default)s)', choices=['PAIRED','SINGLE'])
args=parser.parse_args()

BAM=args.BAM #input MUST be a BWA aligned sorted indexed BAM file
CHR=args.CHR #Y chromosome as defined in the BAM header Y or chrY
len_chr=args.LENGTH #Length of Y chromosome in reference
read_type=args.READ # if single the does not check for proper read pairs use SINGLE
OUT=args.OUTPUT

#PATH TO TOOLS and REFERENCES
#GFile has gene family definition
##START	END	TYPE
##START,END - GENE location based on NCBI RefSeq, or BLAT search for pseudogenes. Make sure there is little overlap between genes
##TYPE - FAMILY NAME or "CONTROL" All the genes within a family should have same name and Control genes should always be named CONTROL(CASE SENSITIVE)
GFile=args.GENE_DEF


#SFile has Y chromosome specific definition HG38 version. (Position	GCpercentage	Mappability	InformativeSites)
##Position - Each nucleotide position on the chromosome of interest(Y)
##GC content - 250bp window with this position as center the percentage of G and C in it.
##Mappability - Mappability score of that position based on 100 base pairs downstream
##Informative Sites - Information of sites that unique to gene family. 
SFile=args.ANNOTATION


#CHECK IF FILES EXIST
if os.path.exists(BAM):
	print "Using "+str(BAM)+" as input BAM file\n"
else:
	print "ERROR: Cannot find input BAM file."
	sys.exit(0)

if os.path.exists(GFile):
	print "Specified path to Gene definition file : "+str(GFile)+"\n"
else:
	print"ERROR: Cannot find Gene definition file "+str(GFile)+" . Please check the code to update the correct path/file name."
	sys.exit(0)

if os.path.exists(SFile):
	print "Specified path to Chromosome Summary file : "+str(SFile)+"\n"
else:
	print"ERROR: Cannot find the file with Chromosome summary file: "+str(SFile)+" with mappability, GCcontent and repeats information. Please check the code to update the correct path/file name."
	sys.exit(0)



def Get_Read_Length(bam_file):
	"""Reads first 1000 reads and obtain the common read length of the sample """
	readlength=[]
	with pysam.AlignmentFile(bam_file, "rb") as bamfile:
		size=1000
		count=0
		for read in bamfile.fetch():
			readlength+=[read.query_length]
			count+=1
			if count>=100:
				break
	return max(set(readlength), key=readlength.count)
			

if Get_Read_Length(BAM) < 100:
	print 'ERROR: Read length of the Input BAM is less than required length of 100 bases. (Tested first 1000 reads)'
	sys.exit(0)

	
#Open the sam file with proper paired reads and filter the reads by alignment
def Count_Matches_CIGAR(cigar_char,cigar_val):
	"""Function to read the parsed CIGAR characters to calculate the number of matches in the first 90 positions """
	i=0 #iterator for while loop
	position=0 #positions in the alignment upto which the cigar string was read
	M=0 #number of matches
	last_val=[] #list to store the position and number of matches before updating in current iteration
	while position < 90 and i <len(cigar_char):
		last_val=[position,M]
		position+=int(cigar_val[i])
		if cigar_char[i]=="M":
			M+=int(cigar_val[i])
		if position>=90: #when the position in alignment crosses the 90 point we want to trim to remove post90 alignment
			extra=position-90 #number of extra positions after 90 
			diff=int(cigar_val[i])-extra #subtract the observed cigar with the extra positions 
			#position=last_val[0]+diff #update the position so we have alignment for first 90 positions
			if cigar_char[i]=="M":
				M=last_val[1]+diff
		i+=1
	return [M,position]

def Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val):
	"""Function to read the parsed MDZtag characters to calculate number of perfect matches in the first 90 positions """
	i=0 #iterator for while loop
	position=0 #positions in the alignment upto which the MDZ tag was read
	MM=0 #number of mismatches and deletions
	while position < 90 and i <len(mismatch_char):
		last_val=[MM,position] #the number of mismatch&deletion before updating in this iteration
		if len(mismatch_char[i])>0:
			if "^" in mismatch_char[i]:
				len_nonM=int(len(mismatch_char[i])-1) #deletion leads with ^, we subtract 1 to ignore the ^
			else:
				len_nonM=int(len(mismatch_char[i]))
		else:
			len_nonM=0
		position+=int(mismatch_val[i])+len_nonM 
		MM+=len_nonM
		if position>=90:
			position=last_val[1]+int(mismatch_val[i])
			if position < 90:
				position=position+len_nonM
				extra=position-90
				diff=len_nonM-extra
				MM=last_val[0]+diff
			else:
				MM=last_val[0]
		i+=1
	return [MM,position]

print "\rFiltering reads for perfect matches"	
#Read the input bam file and parse for proper read pairs and then look for reads with atleast 88 perfect matches in the first 90 base pairs of the read
filtered_read_start_position=BAM+"_"+CHR+"_alignmentSTARTPosition.tab"
with pysam.AlignmentFile(BAM, "rb") as bamfile, open(filtered_read_start_position, "w") as w:
	j,i=0,0
	for read in bamfile.fetch(CHR):
		if read_type=="PAIRED":
			if read.flag == 99 or read.flag == 163 or read.flag == 83 or read.flag == 147:
				cigar_char=re.split('\d+',read.cigarstring)[1:] #parse the character
				cigar_val=re.split('\D+',read.cigarstring)[:-1] #parse the number
				mismatch=read.get_tag('MD')
				mismatch_char=re.split('\d+',mismatch)[1:] #parse the character
				mismatch_val=re.split('\D+',mismatch)	   #parse the number
				if int(cigar_val[0])>=90 and int(mismatch_val[0])>=90:
					w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
				elif int(cigar_val[0])>=90 and int(mismatch_val[0])< 90:
					MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
					if MM[0] <=2 and MM[1]>=90:
						w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")				
				elif int(cigar_val[0])<90:
					M=Count_Matches_CIGAR(cigar_char,cigar_val)
					if M[0]>=88 and M[1]>=90:
						if int(mismatch_val[0])>=88:
							w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
						else:
							MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
							if MM[0] <=2 and MM[1]>=90:
								w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")		
				i+=1
				if i==1000000:
					i=0
					j+=1
					print "\rProcessed "+str(j)+"000000 lines"
		elif read_type=="SINGLE":
			cigar_char=re.split('\d+',read.cigarstring)[1:] #parse the character
			cigar_val=re.split('\D+',read.cigarstring)[:-1] #parse the number
			mismatch=read.get_tag('MD')
			mismatch_char=re.split('\d+',mismatch)[1:] #parse the character
			mismatch_val=re.split('\D+',mismatch)	   #parse the number
			if int(cigar_val[0])>=90 and int(mismatch_val[0])>=90:
				w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
			elif int(cigar_val[0])>=90 and int(mismatch_val[0])< 90:
				MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
				if MM[0] <=2 and MM[1]>=90:
					w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")				
			elif int(cigar_val[0])<90:
				M=Count_Matches_CIGAR(cigar_char,cigar_val)
				if M[0]>=88 and M[1]>=90:
					if int(mismatch_val[0])>=88:
						w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")
					else:
						MM=Count_MisMatch_Deletions_MDZtag(mismatch_char,mismatch_val)
						if MM[0] <=2 and MM[1]>=90:
							w.write(str(read.query_name)+"\t"+str(read.reference_start+1)+"\n")		
			i+=1
			if i==1000000:
				i=0
				j+=1
				print "\rProcessed "+str(j)+"000000 lines"

print "\rFinished filtering reads"
print "\rObtaining chromosome wide read start counts."

#Get the read start counts from alignment start positions
counts={}
with open(filtered_read_start_position,"r") as r:
	for line in r :
		col=line.rstrip('\n').split("\t")
		if col[1] in counts:
			counts[col[1]]+=1
		else:
			counts[col[1]]=1
	


#For every position on Y chromosome, if no reads starting then fill with zeros
#Output file is a list of one-based position specific counts

temp=0
val=[]
for pos in sorted(counts.keys(),key=int):
	diff=int(pos)-temp
	if diff==1:
		val+=[counts[pos]]
	else:
		val+=[0]*(diff-1)
		val+=[counts[pos]]
	temp=int(pos)

#####len_chr=57227415	#Y=57227415
tail_cov=len_chr-temp
val+=[0]*(tail_cov)

####Uncomment the below block and run it if you are trying to debug the later half of the code.
####output file will have the read start counts which were generated from the BAM file.
# output=BAM+"_"+CHR+"_ReadStartCount.txt"
# with open(output, "w") as file:  
	# file.write("StartCount\n")
	# for index in range(len(val)):
		# file.write(str(val[index])+ "\n")

##Get the list of genes whose read depth is needed to be calculated		
print "\rObtaining the gene list"

Gene_list={}
Family_list={}
with open(GFile, "rU") as g: 
	header=g.next()	#START	END	TYPE
	for line in g:
		col=line.rstrip('\n').split("\t") 
		col[2]=col[2].upper()
		Gene_list[str(col[2])+"_"+str(col[0])+"_"+str(col[1])]=col
		if col[2] in Family_list:
			Family_list[col[2]].append(str(col[2])+"_"+str(col[0])+"_"+str(col[1]))
		else:
			Family_list[col[2]]=[str(col[2])+"_"+str(col[0])+"_"+str(col[1])]


##Parse repeat region and add mappability values to RSCcounts
print "\rLoading the Read start counts (RSC) and Mappability values"
Data={} # Dictionary to temporarily store values and convert to data-frame later . Saves runtime
with open(SFile, "rU") as s: 
	header=s.next() #'Position\tGCpercentage\tMappability\tInformativeSites\n'
	for row in s:
		col=row.rstrip('\n').split("\t") 
		#Data[int(col[0])]={'Position': int(col[0]), 'Mappability':float(col[2]), 'GC':float(col[1]), 'RSC':int(val[int(col[0])-1]), 'tInformativeSites':int(col[3]) } #Converting 1 based to 0 based by subtracting 1
		Data[int(col[0])]=[int(col[0]), float(col[2]), float(col[1]), int(val[int(col[0])-1]),int(col[3])]


Summary_data= pandas.DataFrame.from_dict(Data, orient="index")				
Summary_data.columns=['Position', 'Mappability', 'GC', 'RSC','Informative']
Control_data=Summary_data.copy()
Control_data=Control_data.loc[Control_data['Mappability']==1]


Summary_data=Summary_data.values
Control_data=Control_data.values
mean_Control=numpy.mean(Control_data[:,3]) 

print "\rPerforming the GC correction"
###CONTROL & GC correction

#Create 100 windows each having all sites on Y with GC% within that window. Example: window: 0-0.99,1-1.99,2-2.99. All position on Y with GC% >=1 and <2 fall in 1-1.99 window.
GCmean=numpy.empty((0, 1))
for i in range(1,100):
	counts_temp=Control_data[((Control_data[:,2]>=i)==(Control_data[:,2]<(i+1))).nonzero()][:,3]
	if len(counts_temp) > 0 :
		gcm=numpy.mean(counts_temp) #caliculate the mean RSC for each window
	else:
		gcm=0
	#print i, gcm
	GCmean = numpy.append(GCmean,gcm)

#This step below is to make sure there are no 0 to divide with in the next step.
GCmean[GCmean==0]=mean_Control
Correction=(mean_Control/GCmean)
for i in range(len(Correction)):
	if numpy.isnan(Correction[i]):
		Correction[i]=1
	if numpy.isinf(Correction[i]):
		Correction[i]=1

#GC correction step
GCcor_Summary_data=Summary_data.copy()
for i in range(1,100):
	id=((GCcor_Summary_data[:,2]>=(i))==(GCcor_Summary_data[:,2]<(i+1))).nonzero()
	GCcor_Summary_data[id,3]=GCcor_Summary_data[id][:,3]*Correction[i-1]
	

def Control_region_coverage(Y_Summary_data) : # Ychr_summary=['Position', 'Mappability', 'GC', 'RSC','Informative']
		Control_region=Y_Summary_data[(Y_Summary_data[:,1]==1).nonzero()] #All sites with mappability one
		return numpy.mean(Control_region[:,3])
	

def Get_Informative_coverage(Gene_info,Y_Summary_data) : 		# Y_Summary_data=['Position', 'Mappability', 'GC', 'RSC','Informative']; Gene_info=[START,	END,	TYPE]
	Gene_Summary=Y_Summary_data[((Y_Summary_data[:,0]>int(Gene_info[0]))==(Y_Summary_data[:,0]<int(Gene_info[1]))).nonzero()] 		#parse region of the gene (Gene_info[0] is start,Gene_info[1] is end)
	Informative_data=Gene_Summary[(Gene_Summary[:,4]== 1).nonzero()]		#parse informative sites (Informative=1; 0 otherwise)
	return [Gene_info[2],numpy.mean(Informative_data[:,3]),str(Gene_info[2])+"_"+str(Gene_info[0])+"_"+str(Gene_info[1])]  			#return(gene_info,mean coverage)

def Get_ControlGene_coverage(Gene_info,Y_Summary_data) : 		# Y_Summary_data=['Position', 'Mappability', 'GC', 'RSC','Informative']; Gene_info=[START,	END,	TYPE]
	Gene_Summary=Y_Summary_data[((Y_Summary_data[:,0]>int(Gene_info[0]))==(Y_Summary_data[:,0]<int(Gene_info[1]))).nonzero()] 		#parse region of the gene (Gene_info[0] is start,Gene_info[1] is end)
	Unique_data=Gene_Summary[(Gene_Summary[:,1]== 1).nonzero()]		#parse sites with mappability 1
	return [Gene_info[2],numpy.mean(Unique_data[:,3]),str(Gene_info[2])+"_"+str(Gene_info[0])+"_"+str(Gene_info[1])] 

print "\rObtaining the gene level RSC"
Control_coverage=Control_region_coverage(GCcor_Summary_data)

Temp_Coverage={}
for genefamily in Family_list:
	if genefamily == "CONTROL":
		for gene in Family_list[genefamily]:
			Temp_Coverage[gene]=Get_ControlGene_coverage(Gene_list[gene],GCcor_Summary_data)
	elif len(genefamily) == 1 : #If there are single copy genes in gene_def not annotated as CONTROL, but as gene family with one gene
		for gene in Family_list[genefamily]:
			Temp_Coverage[gene]=Get_ControlGene_coverage(Gene_list[gene],GCcor_Summary_data)
	else:
		for gene in Family_list[genefamily]:
			Temp_Coverage[gene]=Get_Informative_coverage(Gene_list[gene],GCcor_Summary_data)


Gene_coverage= pandas.DataFrame.from_dict(Temp_Coverage, orient="index")				
Gene_coverage.columns=['GeneFamily', 'RSCDepth','Gene_ID']

XDG_Genes=Gene_coverage.values[(Gene_coverage.values[:,0]=="CONTROL").nonzero()]
XDG_Control_coverage=numpy.mean(XDG_Genes[:,1])

XDG_CopyNumber=XDG_Genes.copy()
XDG_CopyNumber[:,1]=XDG_CopyNumber[:,1]/Control_coverage

Allgene=Gene_coverage.values
Allgene_CopyNumber=Allgene.copy()
Allgene_CopyNumber[:,1]=Allgene_CopyNumber[:,1]/Control_coverage
Allgene_CopyNumber[Allgene_CopyNumber[:,1].argsort()]

print "\rCaliculating CN values and printing"
Ampliconic_Genes=Gene_coverage.values[(Gene_coverage.values[:,1]!="CONTROL").nonzero()]
AG_CopyNumber=Ampliconic_Genes.copy()
AG_CopyNumber[:,1]=AG_CopyNumber[:,1]/Control_coverage
AG_CopyNumber_XDGbased=Ampliconic_Genes.copy()
AG_CopyNumber_XDGbased[:,1]=AG_CopyNumber_XDGbased[:,1]/XDG_Control_coverage

AG_out=OUT+"Ampliconic_Summary.txt"
with open(AG_out, "w") as AGfile:  
	AGfile.write("GeneFamily\tCopyNumber(MAP=1)\tCopyNumber(XDG)\n")
	for family in sorted(Family_list.keys()):
		if family == "CONTROL": 
			continue
		Family_set=AG_CopyNumber[(AG_CopyNumber[:,0]==family).nonzero()]
		Family_setXDG=AG_CopyNumber_XDGbased[(AG_CopyNumber_XDGbased[:,0]==family).nonzero()]
		#print family+"\t"+str(numpy.sum(Family_set[:,1]))+"\t"+str(numpy.sum(Family_setXDG[:,1]))
		AGfile.write(family+"\t"+str(numpy.sum(Family_set[:,1]))+"\t"+str(numpy.sum(Family_setXDG[:,1]))+ "\n")

AGfile.close()

XDG_out=OUT+"XDG_CopyNumber.txt"
with open(XDG_out, "w") as XDGfile:  
	XDGfile.write("Gene\tCopyNumber\n")		
	for i in XDG_CopyNumber:
		#print i[2]+"\t"+str(i[1])
		XDGfile.write(i[2]+"\t"+str(i[1])+"\n")

# ALL_out=BAM+"AllGenes_CopyNumber.txt"
# with open(ALL_out, "w") as ALLfile:  
	# ALLfile.write("Gene\tCopyNumber\n")		
	# for i in Allgene_CopyNumber[Allgene_CopyNumber[:,1].argsort()]:
		# #print i[2]+"\t"+str(i[1])
		# ALLfile.write(i[2]+"\t"+str(i[1])+"\n")

os.remove(filtered_read_start_position) 