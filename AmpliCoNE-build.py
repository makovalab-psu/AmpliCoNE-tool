import argparse
from Bio import SeqIO
from datetime import datetime as dt
from pprint import pp
import pandas as pd
import subprocess
import pysam
import os
from functools import reduce

class TimeTracker:

    def __init__(self):
        self.start_time = dt.now()
        self.prev_time = dt.now()
    
    def print_time_diff(self, message=None):
        current_time = dt.now()
        overall_time = current_time - self.start_time
        prev_time_diff = current_time - self.prev_time

        print(f"   ----------")
        print(f"   Overall time elapsed: {overall_time}")
        print(f"   Time since last checkpoint: {prev_time_diff}")
        print(f"   Current time: {current_time}")
        if message:
            print(f"   ----------")
            print(f"   message: {message}")
        print(f"   ----------")

        self.prev_time = current_time
    
    def reset(self):
        self.start_time = dt.now()
        self.prev_time = dt.now()


def check_and_get_chr_len(chr_len_file):

    if not os.path.exists(chr_len_file):
        print(f"Error: chromosome length file {chr_len_file} not found")
        exit(1)

    len_chr = 0

    with open(chr_len_file, "r") as f:
        len_chr = int(f.readline().strip())

    if len_chr == 0:
        print(f"Error: chromosome length is 0, check {chr_len_file}")
        exit(1)
    
    return len_chr


def calculate_GC_content(ref, chr_y, tmpdir, outfile):

    if os.path.exists(outfile):
        print(f"found: GC content file {outfile}")
        print("skipping GC content calculation")
        return

    def chunks(seq, window):
        seqlen = len(seq)
        for i in range(0,seqlen,window):
            if (i+window>seqlen):
                j = seqlen 
            else:
                j = i+window
            if i - 250 < 0:
                l = 0
            else:
                l = i - 250
            if j + 250 > seqlen:
                r = seqlen
            else:
                r = j + 250
            yield seq[l:r], i+1, j
            if j == seqlen:
                break

    refFH = open(ref,"r")
    chr_len = 0
    for chr_record in SeqIO.parse(refFH, "fasta"):
        if chr_record.id==chr_y:
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

    if chr_len == 0:
        print(f"Error: chromosome {chr_y} not found in reference: {ref}")
        exit(1)

    with open(outfile, "w") as file:
        file.write("GC\n")
        for index in range(len(GC)):
            file.write(str(GC[index])+ "\n")


def chop_sliding_reference_reads(ref, chr_y, tmpdir, r_len, gene_file, fasta):


    if os.path.exists(fasta):
        print(f"found: fasta file {fasta}")
        print("skipping fasta file generation")
        return


    gene_intervals = []

    with open(gene_file,'r') as f:
        for line in f.readlines():
            line = line.strip()
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            
            if fields[2] == 'CONTROL':
                continue
            
            gene_intervals.append((int(fields[0])-r_len,int(fields[1])+r_len))
            

    # fasta = f"{tmpdir}/Sliding_window{r_len}bp_Reference_reads.fasta"
    #Generate fastq file with each read header has chromosome name and position information (>chr_position) and a sequence is 101bp which starts at position defined in header.
    found = False
    refFH = open(ref,"r")
    for chr_record in SeqIO.parse(refFH, "fasta"):
        print(f"record id {chr_record.id} is it same asi chr_y {chr_y}?")

        if chr_record.id == chr_y:
            
            chr = chr_record.seq
            found = True
            fasta_output=open(fasta,"w")
            for i in range(len(chr)-r_len):
                if i % 1000000 == 0:
                    print(f"Processing position {i} of {len(chr)}")
            
                if any(start <= i <= end for start, end in gene_intervals):
                    r_seq=str(chr[i:i+r_len])
                    r_header = ">" + str(chr_y) + "_" + str(i+1)
                    fasta_output.write(r_header+"\n")
                    fasta_output.write(r_seq+"\n")

            fasta_output.close()

    if not found:
        print(f"Error: chromosome {chr_y} not found in reference: {ref}")
        exit(1)


def parse_informative_sites(SAM, gene_file, debug, r_len, tmpdir, chr_y, informative_sites, informative_list):

    len_chr = check_and_get_chr_len(f"{tmpdir}/{chr_y}_length.txt")

    if os.path.exists(informative_list) and os.path.exists(informative_sites):
        print(f"found: informative sites files {informative_list} and {informative_sites}")
        print("skipping informative sites parsing")
        return
        
    Family_location={}
    with open(gene_file, "r") as g: 

        for line in g:
            if line.startswith("#"):  #START    END    FAMILY
                continue
            col=line.rstrip('\n').split("\t") 
            if col[2] in Family_location:
                Family_location[col[2]].append([col[0],col[1]])
            else:
                Family_location[col[2]]=[[col[0],col[1]]]

    Yposition_mapcount={}
    # print("reading SAM file")
    with pysam.AlignmentFile(SAM, "r") as bamfile:
        for read in bamfile.fetch():
            if read.flag !=4: # 4 = unmapped
                if read.cigarstring == '101M':
                    if read.get_tag("NM")>2:
                        continue
                    if read.query_name in Yposition_mapcount:
                        #print read.query_name, read.flag, read.reference_name, read.reference_start+1, read.cigarstring, read.get_tag("NM")
                        Yposition_mapcount[read.query_name].append(str(read.reference_name)+":"+str(read.reference_start+1))
                    else:
                        Yposition_mapcount[read.query_name]=[(str(read.reference_name)+":"+str(read.reference_start+1))]
            else:
                if read.query_name in Yposition_mapcount:
                    Yposition_mapcount[read.query_name].append('NA')
                else:
                    Yposition_mapcount[read.query_name]=['NA']


    if debug:

        samfile = SAM.split("/")[-1]
        map_file = open(f"{tmpdir}/{samfile}_debug.map", "w")
        pp(Yposition_mapcount, stream=map_file)
        map_file.close()

    #Identify informative sites by parsing reads that mapped to each gene within a family as a set and save it as a dictionary
    # print("Identifying informative sites")
    Family_list_readMapped={}
    for genefamily in Family_location:
        print("---------Gene Family---------")
        print(genefamily)
    
        if genefamily.upper() == 'CONTROL':
            print("Skipping CONTROL gene family")
            continue
        elif len(Family_location[genefamily])==1 : #skip if gene family has only one gene
            print("Skipping gene family with only one gene")
            continue
        else:
            f_size = len(Family_location[genefamily])
            print(f"size: {f_size}")
            print("-----------------------------")
            eachgene_readMapped_perlocation={}
            #print genefamily
            for i in range(f_size):
                gene_start=int(Family_location[genefamily][i][0])
                gene_end=int(Family_location[genefamily][i][1])
                gene_index=range(gene_start,gene_end+1)
                id_gene=str(genefamily)+"_"+str(gene_start)+"_"+str(gene_end)
                #print id_gene
                eachgene_readMapped_perlocation[id_gene]=[]
                temp_list=[]
                for j in gene_index:
                    k=str(chr_y)+"_"+str(j) # NOTE: It should match the read names in SAM file
                    temp_list.append(sorted(Yposition_mapcount[k]))
                eachgene_readMapped_perlocation[id_gene]=set(tuple(i) for i in temp_list)
            Family_list_readMapped[genefamily]=eachgene_readMapped_perlocation

    if debug:

        samfile = SAM.split("/")[-1]
        dict_file = open(f"{tmpdir}/{samfile}_debug.dict", "w")
        pp(Yposition_mapcount, stream=dict_file)
        dict_file.close()

    #For each family label the position that has the number of reads mapped to it equals the size of family and found in all copies as 1 otherwise 0
    Genefamily_output=open(informative_list,"w")
    informative_sites = [0]*len_chr
    for key in Family_list_readMapped:
        #Intersect sets of reads mapped to genes in gene family 
        info_sites=reduce(set.intersection,Family_list_readMapped[key].values())
        size=len(Family_location[key])
        for k in info_sites:
            if len(k)==size:
                x=list(k)
                x.append(key)
                Genefamily_output.write(str('\t'.join(x)+"\n"))
                for l in k:
                    informative_sites[int(l.split(':')[1])-1]=1 #Converting read location 1 based to python 0 based by subtracting one here

    Genefamily_output.close()


    informative_output=open(informative_sites,"w")
    informative_output.write("Informative_Sites\n")
    for index in range(len(informative_sites)):
        informative_output.write(str(informative_sites[index])+ "\n")

    informative_output.close()


def map_fasta_to_sam(fasta, sam, tmpdir):

    # check if bowtie-index is present
    run_bowtie2_index = True
    if os.path.exists(f"{tmpdir}/bowtie2_index.1.bt2"):
        print("found: bowtie2 index in tmp directory")
        print("skipping bowtie2 index creation")
        run_bowtie2_index = False


    if os.path.exists(sam):
        print("found: sam file in tmp directory")
        print("skipping mapping")
        return

    if run_bowtie2_index:
        print("creating bowtie2 index")
        result = subprocess.run(["bowtie2-build", fasta, f"{tmpdir}/bowtie2_index"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        print(result.stdout)
        print(result.stderr)

    print("mapping fasta to sam")
    result = subprocess.run(["bowtie2", "-k", "15", "--threads","32","-x", f"{tmpdir}/bowtie2_index", "-f", fasta, "-S", sam], stdout=subprocess.PIPE, stderr=subprocess.PIPE)


def parse_repeatmasker_output(repeatmasker, chr_y, tmpdir, outfile):


    len_chr = check_and_get_chr_len(f"{tmpdir}/{chr_y}_length.txt")
   
    if os.path.exists(outfile):
        print(f"found: repeat masked file {outfile}")
        print("skipping repeatmasker parsing")
        return


    repeatmask_val = [1]*len_chr
    with open(repeatmasker,"r") as repeatmasker_file:

        for line in repeatmasker_file:
            if line.startswith("#"):
                continue
            if line.strip() == "":
                continue
            col=line.rstrip('\n').split()
            if col[4] == chr_y:
                repeatmask_val[int(col[5]):int(col[6])]=[0]*(int(col[6])-int(col[5]))

    file = open(outfile, "w")
    file.write("RepeatMasked\n")
    for index in range(len(repeatmask_val)):
        file.write(str(repeatmask_val[index])+ "\n")

    file.close()


def parse_tandem_repeat_finder_output(tandem_repeat_finder, chr_y, tmpdir, outfile):

    len_chr = check_and_get_chr_len(f"{tmpdir}/{chr_y}_length.txt")

    if os.path.exists(outfile):
        print(f"found: tandem repeat masked file {outfile}")
        print("skipping tandem repeat finder parsing")
        return

    trf_val = [1]*len_chr
    with open(tandem_repeat_finder,"r") as trFH:
        for line in trFH :
            col=line.rstrip('\n').split('\t')
            if col[0] == chr_y:
                trf_val[int(col[1]):int(col[2])]=[0]*(int(col[2])-int(col[1]))


    trFH.close()

    file = open(outfile, "w")
    file.write("TRF\n")
    for index in range(len(trf_val)):
        file.write(str(trf_val[index])+ "\n")

    file.close()


def parse_gem_mappability_output(gem_mappability, chr_y, tmpdir, mappability_file):

    len_chr = check_and_get_chr_len(f"{tmpdir}/{chr_y}_length.txt")
    outputfile = f"{tmpdir}/{chr_y}_MAPPABILITY_bypos.txt"

    if os.path.exists(outputfile):
        print(f"found: mappability file {outputfile}")
        print("skipping mappability parsing")
        return

    df = pd.read_table(gem_mappability, sep="\t", header=None, low_memory=False)
    chrY=df[(df[0]==chr_y)] #parse Y chromsome specific reads
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
    with open(outputfile, "w") as file:
        file.write("Mappability\n")
        for index in range(len(map_val)):
            file.write(str(map_val[index])+"\n")


def generate_summary_information(GC, MAP, RM, TRF, INFO, OUT, tmpdir, chr_y):

    len_chr = check_and_get_chr_len(f"{tmpdir}/{chr_y}_length.txt")

    if os.path.exists(OUT):
        print(f"found: output file {OUT}")
        print("skipping summary information generation")
        return

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


def main():

    parser=argparse.ArgumentParser(
        description='''AmpliCoNE preprocessing.''')
    parser.add_argument('--chromosome','-c', help="chromosome name as in reference", required=True )
    parser.add_argument('--reference', '-i', help="reference file", required=True) # i for consistency with AmpliCoNE v1.0.0
    parser.add_argument('--gem-mappability','-m', help="GEM-MAPPABILITY output in bed format", required=True)
    parser.add_argument('--repeat-masker','-r', help="RepeatMasker output", required=True)
    parser.add_argument('--tandem-repeat-finder','-t', help="Tandem Repeat Finder output", required=True)
    parser.add_argument('--gene-definition','-g', help="gene definition file", required=True)
    parser.add_argument('--output','-o', help="output file name", required=True)
    parser.add_argument('--tmpdir','-d', help="temporary directory", default="/tmp")
    parser.add_argument('--debug','-D', help="debug mode", action="store_true")
    parser.add_argument('--read-length','-l', help="read length", default=101)
    parser.add_argument('--step','-s', help="run single step of the build process", default=0)

    args = parser.parse_args()

    chr_y = args.chromosome
    reference = args.reference
    gem_mappability = args.gem_mappability
    repeat_masker = args.repeat_masker
    tandem_repeat_finder = args.tandem_repeat_finder
    gene_definition = args.gene_definition
    output = args.output
    tmpdir = args.tmpdir
    debug = args.debug
    r_len = args.read_length
    step = 0
    step = int(args.step) 

    time_tracker = TimeTracker()


    # these files are created in the individual steps
    # if they are present, individual steps can be run
    gc_file = f"{tmpdir}/{chr_y}_GC_Content.txt"
    fasta = f"{tmpdir}/Sliding_window{r_len}bp_Reference_reads.fasta"
    sam = f"{tmpdir}/Sliding_window{r_len}bp_Reference_reads.sam"
    repeat_masked = f"{tmpdir}/{chr_y}_REPEAT_MASKED.txt"
    tandem_repeat_masked = f"{tmpdir}/{chr_y}_TANDAMREPEAT_MASKED.txt"
    informative_sites = f"{tmpdir}/{chr_y}_Informative_{r_len}bp.tab"
    mappability_file = f"{tmpdir}/{chr_y}_MAPPABILITY_bypos.txt"
    informative_list = f"{tmpdir}/{chr_y}_Informative_list_GeneFamily.tab"


    #test if bowtie is installed
    try:
        result = subprocess.run(["bowtie2", "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    except FileNotFoundError:
        print("Error: bowtie2 not found")
        exit(1)


    if not os.path.exists(tmpdir):
        os.makedirs(tmpdir)


    if step == 0 or step == 1:
        if debug:
            time_tracker.print_time_diff("Calculate GC content")
        calculate_GC_content(reference, chr_y, tmpdir, gc_file)

    if step == 0 or step == 2:
        if debug:
            time_tracker.print_time_diff("Create Fasta file with chopped sliding reference reads")
        chop_sliding_reference_reads(reference, chr_y, tmpdir, r_len, gene_definition, fasta)

    if step == 0 or step == 3:
        if not os.path.exists(fasta):
            print("Cannot run step 3 Mapping without fasta file")
            print(f"Error: fasta file {fasta} not found")
            exit(1)
        if debug:
            time_tracker.print_time_diff("Create SAM file by mapping fasta file to reference")
        map_fasta_to_sam(fasta, sam, tmpdir)

    if step == 0 or step == 4:
        if debug:
            time_tracker.print_time_diff("Parse RepeatMasker output")
        parse_repeatmasker_output(repeat_masker, chr_y, tmpdir, repeat_masked)

    if step == 0 or step == 5:
        if debug:
            time_tracker.print_time_diff("Parse Tandem Repeat Finder output")
        parse_tandem_repeat_finder_output(tandem_repeat_finder, chr_y, tmpdir, tandem_repeat_masked)

    if step == 0 or step == 6:
        if debug:
            time_tracker.print_time_diff("Parse GEM-Mappability output")
        parse_gem_mappability_output(gem_mappability, chr_y, tmpdir, mappability_file)
    
    if step == 0 or step == 7:
        if debug:
            time_tracker.print_time_diff("Parse informative sites")
        parse_informative_sites(sam, gene_definition, debug, r_len, tmpdir, chr_y, informative_sites, informative_list)

    if step == 0 or step == 8:
        if debug:
            time_tracker.print_time_diff("Create summary information")
        generate_summary_information(gc_file, mappability_file, repeat_masked, tandem_repeat_masked, informative_sites, output, tmpdir, chr_y)


    print("Done")

if __name__ == "__main__":
    main()