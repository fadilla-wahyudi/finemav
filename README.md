# *FineMAV*
This stand-alone program implements the Fine-Mapping of Adaptive Variation (*FineMAV*) statistic for detection of positively selected variants.

The *FineMAV* score of the derived allele for each SNP is calculated by multiplying three metrics:

- Derived allele purity -- measure of population differentiation that ranges from 0 to 1, where a score of 1 indicates that all the derived alelles fall into a single population. 
- Derived allele frequency 
- The Combined Annotation-Dependent Depletion PHRED-scaled C-score (CADD_PHRED) -- a measure of functionality

-**Add screenshot of genome-wide scores on IGV**
-**Zarul said add explanation on how to use genome browsers... perhaps consider?**

## Table of Contents
- [*FineMAV*](#finemav)
- [Installation](#installation)
  - [GitHub](#github)
  - [Anaconda](#anaconda)
  - [GUI version](#gui-version)
- [Running *FineMAV*](#running-finemav)
  - [Checklist of fields that need to be provided](#checklist-of-fields-that-need-to-be-provided)
  - [Usage](#usage)
  - [List of required and optional flags](#list-of-required-and-optional-flags)
- [Extracting information that is not annotated in the VCF file](#extracting-information-that-is-not-annotated-in-the-vcf-file)
  - [AF](#af)
  - [AA and CADD_PHRED](#aa-and-cadd_phred)
- [Default penalty parameters](#default-penalty-parameters)
- [Citation](#citation)
- [License](#license)
  
## Installation

### GitHub
First, clone this repository.
```
git clone https://github.com/fadilla-wahyudi/finemav/
```
Give the binary file executable permissions.
```
cd finemav
chmod +x finemav
```
To access the file globally, move it to one of the directories in your PATH. To check the PATH variable, type in `echo $PATH`.

### Anaconda
- **Put it on anaconda when the paper goes through**

### GUI version
A GUI version can be downloaded here: https://drive.google.com/file/d/1xBhQGpUhVd02kyIuevVIuqac4zJ_13Tm/view?usp=sharing

- **Change the link when published**
- **Need to edit the the appearance**
- **Add screenshot**


## Running *FineMAV*

### Checklist of fields that need to be provided
For the program to work, the input file (e.g. *VCF_AllFields.txt*) should contain the following fields in a tab-delimited format with the field names as the header:

Fields | Description | Mandatory VCF column
-------|-------------|--------------------
CHROM:POS | Chromosome number:position | Yes
ID | Identifier | Yes
REF | Reference allele | Yes
ALT | Alternative allele (this version can only use biallelic SNPs) | Yes
AA | Ancestral allele | No
CADD_PHRED | PHRED-scaled Combined Annotation Dependent Depletion (CADD) score | No
AF | Allele frequency for the ALT allele (this is done for each population) | No

In an ideal situation, all these fields can be found in the VCF file in which the non-mandatory VCF columns are annotated in the INFO field. If that is the case, the information can be extracted using software like BCFtools:

```
bcftools query --format '%CHROM:%POS\t%ID\t%REF\t%ALT\t%INFO/AA\t%INFO/CADD_PHRED\t%INFO/AF\n' --print-header input.vcf.gz > output.txt
```
If these fields are not found in the VCF file, it can be [retrieved by other means](#extracting-information-that-is-not-annotated-in-the-vcf-file). 

### Usage
We have created test files for you to try out. Here's an example on how to use one of them:

```
finemav --input-file test_files/VCF_AllFields.txt --prefix testrun --reference-genome hg19
```
The input file (*VCF_AllFields.txt*) consists of all the fields mentioned above in a tab-delimited format. It comprises of 11 SNPs and three populations named POP1, POP2 and POP3.

There should be three types of output files generated:
1. Log file (*testrun_finemav.log*)
2. bigWig files that can be loaded onto a genome browser for visualisation (*testrun_POP1.bw*, *testrun_POP2.bw* and *testrun_POP3.bw*)
3. Tab-delimited table containing the *FineMAV* scores and the intermediate calculations (*testrun.txt*)

The columns in the tab-delimited output table (*testrun.txt*) is as follows:

Fields | Description 
-------|-------------
CHROM | Chromosome number
POS | Position
ID | Identifier
REF | Reference allele
ALT | Alternative allele
AA | Ancestral allele
CADD_PHRED | PHRED-scaled CADD score
AF_*** | Allele frequency for POP1, POP2 and POP3
DER | Derived allele
DAF_*** | Derived allele frequency for POP1, POP2 and POP3
SUM_OF_DER | Sum of the derived allele frequency of all populations
DAP | Derived allele purity
FineMAV_*** | *FineMAV* scores for POP1, POP2 and POP3


### List of required and optional flags
Required flags
- `-i --input-file` // this is the tab-delimited file containing the fields that can be extracted from the VCF file
- `-x --prefix` // prefix assigned to the output files
- `-r --reference-genome` // GRCh37/hg19 or GRCh38/hg38

Optional flags
- `-v --vep-file` // if AA and/or CADD_PHRED is not annotated in the VCF file, there is an option to supplment that information using Ensembl's Variant Effect Predictor (VEP). To find out how more on how to do that, go to the [AA and CADD_PHRED](#aa-and-cadd_phred) section. 
- `-p --penalty` // specify the penalty parameter if you do not want to use the [default parameters](#default-penalty-parameters)
- `-c --chunksize` // specifies the number of lines per chunk (default=200000)

## Extracting information that is not annotated in the VCF file

### AF
If the AF is not annotated in the VCF file, we suggest using the BCFtools [+fill-tags plugin](https://samtools.github.io/bcftools/howtos/plugin.fill-tags.html) where it can calculate the AF for each population. The +fill-tags plugin can compute the AF in the INFO column of the VCF file when it is provided with the list of samples and which population they belong to. The command can be piped to the [BCFtools "query" command](http://samtools.github.io/bcftools/bcftools.html#query) so that it can be extracted together into a tab-delimited file.

Here is an example of a BCFtools computing the AF of two populations, POP1 and POP2:
```
bcftools +fill-tags input.vcf.gz -- --tags AF --samples-file list_of_samples.txt | bcftools query --format '%CHROM:%POS\t%ID\t%REF\t%ALT\t%INFO/AA\t%INFO/CADD_PHRED\t%INFO/AF_POP1\t%INFO/AF_POP2\n' --print-header > output.txt

```

### AA and CADD_PHRED
This tool allows users to supplement the AA and CADD_PHRED information using the `--vep-file` flag if these fields are missing from the VCF file. AA and CADD_PHRED should be annotated using the command line [Ensembl VEP](https://asia.ensembl.org/info/docs/tools/vep/script/index.html)'s [AncestralAllele plugin](https://github.com/Ensembl/VEP_plugins/blob/release/101/AncestralAllele.pm) and [CADD plugin](https://github.com/Ensembl/VEP_plugins/blob/release/101/CADD.pm), respectively. The `--vep-file` should be a tab-delimited file.

Here's an example of a VEP command used to annotate a VCF file's AA and CADD_PHRED information: 
```
vep -i input.vcf.gz --plugin CADD,whole_genome_SNVs.tsv.gz --plugin AncestralAllele,homo_sapiens_ancestor_GRCh37_e71.fa.gz --tab --fields "Location,Uploaded_variation,AA,CADD_PHRED" -o output.txt
```
It is not required for the "Uploaded_variation" column, which indicates the SNP ID, to be included in the VEP-generated file. You can chose to include the "Uploaded_variation" if you feel that the "ID" column in the VCF file is not updated. If you chose to include the "Uploaded_variation" column, the "ID" field should be excluded from the `--input-file` file. 

We have created test files for you to try out. These test files take into account different scenarios where the AA and/or CADD_PHRED fields may be missing from the VCF file.

Supplemented using VEP | Command
-----------------------|--------
AA | `finemav -i VCF_no_AA.txt -x output -r hg19 -v VEP_AA.txt`
CADD_PHRED | `finemav -i VCF_no_CADDPHRED.txt -x output -r hg19 -v VEP_CADDPHRED.txt`
AA, CADD_PHRED | `finemav -i VCF_no_AA_CADDPHRED.txt -x output -r hg19 -v VEP_AA_CADDPHRED.txt`


## Default penalty parameters
The derived allele purity metric relies on a penalty parameter to maximise the magnitude between positively selected and nearby neutral variants.
The suggested default penalty parameters seen below were determined empirically by [Szpak et al.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1380-2). 


Number of populations | Penalty
----------------------|--------
2 | 4.962
3 | 3.500
4 | 2.981
5 | 2.271
6 | 2.533
7 | 2.411

## Citation
If you utilise this program, please cite the following papers together:
>Szpak, M., Mezzavilla, M., Ayub, Q. et al. *FineMAV*: prioritizing candidate genetic variants driving local adaptations in human populations. Genome Biol **19**, 5 (2018). https://doi.org/10.1186/s13059-017-1380-2

## License
- **Create MIT License**
- **Read more about licensing a GitHub repository here: https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/licensing-a-repository**
