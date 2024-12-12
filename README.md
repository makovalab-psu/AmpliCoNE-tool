# AmpliCoNE: Ampliconic Copy Number Estimator

A tool to estimate the copy number of ampliconic gene families in human Y chromosome using Illumina whole genome sequencing data.


## Installation

### Dependencies
For the latest version of AmpliCoNE we recommend using Python 3.9 or later. The tool was tested and works well with the following list of tools and packages:
```
Python 3.12.x
Numpy	2.1.2
Pandas 2.2.3	
Pysam	0.22.1
Biopython 1.84
Bowtie2 version 2.4.2 
```

#### Steps to install the dependencies using conda 

Create an environment.yml file as described below.
```
vim amplicone-environment.yml

name: amplicone
channels:
  - bioconda
dependencies:
  - python=3.12
  - pandas=2.2.3
  - numpy=2.1.2
  - pysam=0.22.1
  - biopython=1.84
  - bowtie2=2.*
```
Use the amplicone-environment.yml file to create an environment named "amplicone"
```
conda env create -f amplicone-environment.yml
```
Load the environment each time you run AmpliCoNE.
```
source activate amplicone
```

### Install AmpliCoNE

Once the above mentioned dependencies are installed, clone the AmpliCone tool from the repository https://github.com/makovalab-psu/AmpliCoNE-tool by running:
```
git clone https://github.com/makovalab-psu/AmpliCoNE-tool.git
```



## AmpliCoNE usage for human chromosome Y ampliconic genes using hg38

### Download the Y chromosome annotation file and the gene definition file [here](http://www.bx.psu.edu/medvedev_lab/amplicone/hg38/hg38_amplicone_files.tar.gz). 

```
#download annotation and gene definition files for hg38
wget http://www.bx.psu.edu/medvedev_lab/amplicone/hg38/hg38_amplicone_files.tar.gz

#uncompress the files
tar -zxvf hg38_amplicone_files.tar.gz

#list files downloaded (gene_definition_hg38.tab  hg38_Ychromosome_annotation.tab)
ls hg38_amplicone_files/

```

### Generate input BAM files

Generate BAM files by aligning your Illumina whole genome sequencing dataset to hg38 using BWA-mem. (hg38 version reference can be downloaded from [UCSC genome browser](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/) or [ENSEMBL](http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/)). The BAM files must be sorted by position and indexed.


Download test dataset [here](http://www.bx.psu.edu/medvedev_lab/amplicone/hg38/test.tar.gz).
```
#download test BAM file
wget http://www.bx.psu.edu/medvedev_lab/amplicone/hg38/test_data.tar.gz

#uncompress the files
tar -zxvf test_data.tar.gz

#list files downloaded (test.bam  test.bam.bai)
ls test_data/

```

### AmpliCoNE-count usage

AmpliCoNE-count.py takes the gene definition file, annotation file and BAM file as input to estimate the copy number of the nine ampliconic gene families. In addition, the name of the Y chromosome as defined in reference and its sequence length are required parameters. 

```
python AmpliCoNE-count.py --GENE_DEF hg38_amplicone_files/gene_definition_hg38.tab --ANNOTATION hg38_amplicone_files/hg38_Ychromosome_annotation.tab --BAM test_data/test.bam --CHR Y --LENGTH 57227415
```

#### Other parameters in AmpliCoNE-count.py
The length of the chromosome must be set while using AmpliCoNE for a reference other than hg38. 
```
--LENGTH <int> (default:57227415; length of Y chromosome in hg38) 
```
IMPORTANT: DO NOT USE DEFAULT VALUE FOR NON HG38 REFERENCE

Parameter to define if the BAM file contains single end reads.

```
--READ PAIRED or SINGLE (default:PAIRED)
```

## Output description

AmpliCoNE-count will generate two output files:
 - \<OUTPUT>Ampliconic_Summary.txt
  
    A tab separated file with the ampliconic gene family copy number. First column will have the family name. Second column will have       the gene copy number estimated using all the sites on Y with mappability 1. Third column will have the gene copy number estimated     using the X-degenerate genes (CONTROL) defined in gene definition file.
    
    For example:
 
    | GeneFamily | CopyNumber(MAP=1) | CopyNumber(XDG) |
    |----------- |:-----------------:|----------------:|
    | BPY2       |       3           |       3         |
    | CDY        |       4           |       4         |
    | DAZ        |       4           |       4         |
    | HSFY       |       2           |       2         |
    | PRY        |       2           |       2         |
    | RBMY       |       7           |       7         |
    | TSPY       |      22           |      22         |
    | VCY        |       4           |       4         |
    | XKRY       |       2           |       2         |
    
  
  - \<OUTPUT>XDG_CopyNumber.txt
  
    A tab separated file with two columns. First column will have the X-degenerate gene (XDG) ids. Second column will have the gene copy     number estimates.
    
    For example:
    
    | Gene       | CopyNumber(MAP=1) | 
    |----------- |:-----------------:|
    | GeneID1    |       1           |
    | GeneID2    |       1           |
    | GeneID3    |       1           |
    | GeneID4    |       1           |
    
    XDG_CopyNumber.txt file can provide information about the quality of sample. We can check if the estimates of the CONTROL genes copy number is close to one. 



## AmpliCoNE usage with other reference genomes / species

Before running the tool on a dataset from a different species, you must first perform steps 1-3 below. 
This needs to be run only once, and the resulting files can be reused when analyzing datasets from the same species. 
After these steps, AmpliCoNE-count can be run for each sample to estimate copy number, as described above.


### Step 1: Download pre-requisite files

The pre-requisite files for most genomes can be downloaded from UCSC genome browser. Make sure that the chromosome annotation (chrY, Y) is the same in all the files. If the files are not available for download, please generate them using  GEM-library, RepeatMasker, and TandemRepeatFinder tools. 


- Reference genome (Example : hg38, FASTA format)

  Note: The Y chromosome should be present as one continuous scaffold in the reference file. 

- Reference specific mappability file (Generated using GEM library; BED format)
```
#Steps to generate mappability
gem-indexer -i <REF.fa> -o <REF.fa> --complement emulate --verbose
gem-mappability -I <REF.fa.gem> -l 101 -o <OUT_MAPPABILITY> -m 2 -e 2
gem-2-wig -I <REF.fa.gem> -i <OUT_MAPPABILITY.mappability> -o <WIG_OUT_MAPPABILITY>
wig2bed <WIG_OUT_MAPPABILITY.wig> <REF_MAPPABILITY.bed>
```

- Reference specific RepeatMasker output (.out file)

- Reference specific Tandem Repeat Finder output (BED format)

Note: RepeatMasker and Tandem Repeat Finder only need to be run on the Y chromosome. 

### Step2: Generate gene definition file


- Format  : TSV
- Columns : Start Position, End Position and Gene Family Name
- Rows    : Gene copy

NOTE: For gene copies in control, the gene family name (3rd column) should be "CONTROL".(case sensitive)

**How to generate the file**:
- Identify a representative gene for each gene family and extract its sequence from NCBI. (E.g. most common variant or consensus sequence built from individual gene copies within the gene family.)
- BLAT the representative gene sequence against reference genome to identify all locations with greater than 99% identity.
- Represent each location identified with its start and end location followed by the gene family name.

NOTE: All the locations >99% identity representing a gene family should be present on Y chromosome only. Each gene copy within a family should be represented as a row in the gene definition file. 

For human ampliconic genes the gene definition file is available [here](http://www.bx.psu.edu/medvedev_lab/amplicone/hg38/hg38_amplicone_files.tar.gz). Below is an example file:

| START   |  END   | TYPE    |
|-------- |:------:|--------:|
| 100     |  300   | GF1     |
| 750     |  1050  | GF1     |
| 3000    |  3300  | GF1     |
| 9900    |  10900 | GF2     |
| 2220    |  2320  | GF2     |
| 3500    |  3750  | CONTROL |
| 4190    |  4420  | CONTROL |
| 4770    |  5040  | CONTROL |


### Step 3: Generate Ychr_annotation file using AmpliCoNE-build

```
sh AmpliCoNE-build.py -c <chrY> -i <REF.fasta> -m <REF_MAPPABILITY.bed> -r <REF_REPMASK.out> -t <REF_TRF.bed> -g <Gene_Definition.tab> -o chrY_annotation.tab
```

Use the output file (annotation file) from AmpliCoNE-build.sh and gene definition file to run AmpliCoNE-count.py and estimate the gene family copy numbers.

```
python AmpliCoNE-count.py --GENE_DEF <Gene_Definition.tab> --ANNOTATION <chrY_annotation.tab> --BAM <BAM> --CHR <chr> --LENGTH <int>
```
