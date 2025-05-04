# TE-discovery
## Description
These scripts are a pipeline used to identify novel transposable element sequences from eukaryotic genomes that have the structure:

![](images/Transposon_picture.png)

The innovation of these scripts compared to alternative methods (e.g. RepeatModeler) is that they identify possible transposable elements using structure rather than sequence homology (allowing brand new elements to be discovered), and that manual (human) review steps are incorporated into the pipeline to ensure high quality results. These scripts emphasize quality over quantity. Therefore, some elements may not be detected, but those that are detected (even if there only a few copies in the genome) are unlikely to be artifacts. 

The pipeline is run as a series of steps. Each can step can be run individually or in sequence, depending on the computer resources available for each step.

***STEP 1 Map proteins to the genome*** 

Map all the protein sequences to the genome and filter out protein sequences that either: 1) have low copy numbers (you'd expect a transposable element sequence to be repeated), 2) don't map with sufficient length and identify, 3) only be made up of a single cluster of locally duplicated copies. These steps are inspired from the transposable element discovery pipeline described by [Goubert et al. 2022](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-021-00259-7).

***STEP 2 Scan around the proteins for TSD-TIR patterns***
Scan around each protein for properly located TSDs and TIRs. All the genomic copies are kept as as a group and aligned. 

***STEP 3 Manual review of potential transposable elements***
The TSD-TIR junctions of each copy of each potential element are displayed in human readable format to verify that the structure appears to be a genuine transposable element.

## Running the scripts
### Input Sequences
At its core any analysis requires two files 1) a genome sequence and 2) a list of protein sequences within that genome that are potentially those transposable elements. The second file could be, for example, a list of all identified protein sequences from the genome.

### Running the pipeline 

**Required of all steps**

The pipeline is written in Perl, and is executed using 

	perl mainscript.pl -n <folder name> -e <folder name> <other parameters>

Where:

	-n folder name where to store (or where have already been stored) analysis files, these are input parameter lists and reports on the analysis
	-e folder name where to store (or where have already been stored) transposable elements sequences and related information

and other parameters are step specific.

**Specify steps to execute**

Unless otherwise specified all the analysis steps above will be executed. To specify that just one or a range
of consecutive steps be executed set the parameter(s):

    -a number of first (or only) step to execute
    -b number of last step to execute, if not specified then only the first step will be executed

### Step specific parameters

***STEP1 parameters*** 

The input file(s) for this step can be either a fasta file of proteins and an associated genome file, or a tblastn output file.
The tblastn option is so that this time consuming step can be run on a faster machine on its own (e.g. a computer cluster), the tblastn output must be formatted in the same way as in the script. Specifically, the `tblastn` output parameter must be set to:

	-outfmt "6 qseqid sseqid sstart send pident length qlen"

**If STEP1 uses a list of proteins as input** specify the following parameters:
 
    -p fasta file of proteins that could be part of a TE 
    -g fasta file of the genome associated with the proteins; this file must have been formatted with "makeblastdb" ahead of time

**If STEP1 uses a tBLASTN file as input** specify the following parameter:
 
    -t tblastn output in the same format used in the script 

***STEP2 parameter and file requirement*** 

When running this step the `-g` parameter as described in STEP1 must be specified. In addition, the script will expect a file to be present in the same directory as the genome file. This file must have the name of the genome file followed by `.length`. Each line in this file should contain a fasta title from the genome file and the length in nucleotides of the associated sequence. This `.length` file can be generated using the following BASH commands:

	samtools faidx <genome fasta formatted file name>
	awk '{OFS="\t"; print \$1,\$2}\' < <genome fasta formatted file name>.fai > <genome fasta formatted file name>.length
    
***STEP3 has no step specific parameters or requirements*** 

## Understanding analysis and README files

Each step will generate files that specify the outcome of the analysis. The goal of these files is so that the fate of every input sequence can be traced and, if applicable, the exact reason why it was rejected can be determined.

### README files
Every protein sequence that was considered as a transposable element (i.e. it was not filtered out during STEP1) has its own folder. These folders are found in one of two locations: 1) inside the folder specified by the `-e` parameter, these are the best candidates for real transposable elements; 2) inside the folder specified by the `-n` parameter, in the sub-folder `Rejected_elements`, these have been excluded from consideration. Inside each of these sub-folders is a `README.txt` file that may contain the following information (depending how long it was considered a candidate for being a transposable element).

1) The number of genomic locations matching this sequences using tBLASTN
2) Name of .bed file with blast locations
3) Name of .alipos file that contains information on element alignment locations removed because they are shared by too few copies of this element (this is the same proceedure as described in [Goubert et al. 2022](https://mobilednajournal.biomedcentral.com/articles/10.1186/s13100-021-00259-7). This file may be used to reconstruct the original alignment, prior to trimming the low copy position.
4) Name of .maf file that contains the multiple sequence alignment of all the extended blast matches (see STEP2). Alignment positions with too many gaps are removed from this alignment (see .alipos file).  
5) A statement of how the element performed on the "edge test". A true transposable element is expected to have a sharp transition in the alignment between sequences outside the element and TSD-TIR sequences inside the element. If such a sharp transition is too close the edge of the extended element, the script assumes no such transition is present and the sequences are not considered to be from transposable elements.
6) Name of .tirtsd file that is meant to be used with the .maf file. This .tirtsd file specifies alignment positions that delimit possible outside TIR positions (in the .tirtsd file these are `loc1` and `loc2`). Each line also specifies a success code (see below), the number of alignment sequence that have intact TIRs, and the number of intact TSDs (TSD categories are TA, 2bp., 3bp, ... , 10bp.). The README also specifies the number of alignment sequences that were rejected because their success code excluded them from being considered as transposable elements.
7) A statement from the manual review in STEP3, including notes from the reviewer.

print ANALYSIS "STEP2: TIR-TSD success code, first number: Is the proportion of sequences with TIRs higher than MIN_PROPORTION_SEQ_WITH_TIR?\n";
print ANALYSIS "STEP2: TIR-TSD sode, second number: Does this candidate have one of the highest number of TIRs?\n";
print ANALYSIS "STEP2: TIR-TSD sode, third number: Do the TIRs tend to start and end with the same bases?\n";
print ANALYSIS "STEP2: TIR-TSD sode, first number: Does this candidate have one of the highest number of TSDs?\n";
    

