# Feb 2025 - Top-level script for TE discovery pipeline
# Inputs: List of TE sequences (FASTA) and a genome (FASTA)

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use FindBin::libs;  # this sets up that the directory lib will have all modules necessary to run the program 
use TEdiscovery;    # these are the subcripts necessary to run the pipeline

### CONSTANTS
my $NUM_THREADS = 8; # number of threads to use

### INPUTs from command line
my $INPUT_PROTEIN_SEQUENCES; # fasta formated file with input protein sequences
my $INPUT_GENOME; # fasta formated file with genome that input proteins
my $ANALYSIS_NAME; # name to give to this analysis, can be any string
my $TEMP_OUTPUT_DIR; # directory where analsysis output files of the analysis are stored, these are not destined to be permananet
my $START_STEP = 0; # analysis step to start at, by default is set to zero
my $END_STEP = 100000; # analysis step to end at, by deault set to 1000000 (hopefully fewer than those number of steps)

### CHECK INPUTS Read and check that the inputs have been provided
GetOptions(
	't:s'   => \$INPUT_PROTEIN_SEQUENCES,
	'g:s'   => \$INPUT_GENOME,
    'n:s'   => \$ANALYSIS_NAME,
    'a:i'   => \$START_STEP,
    'b:i'   => \$END_STEP
);

## CHECK INPUTS Validate required input files
unless (\$INPUT_PROTEIN_SEQUENCES and $INPUT_GENOME and $ANALYSIS_NAME) {
	die "usage perl TEdiscovery.pl <-t fasta file with starting TE sequences REQUIRED> <-g fasta file of genome processed with makeblastdb under the same name REQUIRED> <-n name for this analysis REQUIRED> <-a start the analysis at this step OPTIONAL> <-b end the analysis at this step OPTIONAL\n";
}

## CHECK INPUTS Create temporary output directory if necessary
my $TEMP_OUTPUT_DIR="./$ANALYSIS_NAME-outfiles"; # directory to store output files of current analysis (no slash at the end), these can be destroyed when analysis is finished
if (-d $TEMP_OUTPUT_DIR) {
    warn "Directory $TEMP_OUTPUT_DIR already exists, using it for this analysis\n";
}
else {
    warn "Creating directory $TEMP_OUTPUT_DIR for temporary storage\n";
    mkdir( $TEMP_OUTPUT_DIR ) or die "Couldn't create $TEMP_OUTPUT_DIR directory, $!";
}

### PIPELINE STEP 1 identify proteins that match the genome in multiple locations
my $step_number = 1;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## Excute the tblastn search
    `tblastn -query $INPUT_PROTEIN_SEQUENCES -db $INPUT_GENOME -outfmt "6 qseqid sseqid sstart send qseqid pident length qlen" -out $TEMP_OUTPUT_DIR/tblastn.o -num_threads $NUM_THREADS`;
    if ($?) { die "ERROR executing tblastn, stoping analysis: error code $?\n"}

    ## Inspired by the Goubert et al. protocol, filter elements that 1) have >= 80% identity to genome, 2) have 50% length of the query, 3) occur multiple times
    my $genome_identity = 80; # per protein, minimum percent identity between protein and genome
    my $coverage_ratio = 0.5; # per protein, minimum ratio of (blast match length) / (query length)
    my $copy_number = 2; # minimum number of copies that hit different parts of the genome 

    ### need keep only one copy of overlaping hits
    
#    `#awk '{OFS="\t"; if ($3 >= 80 && (($4/$13) > 0.5 )) {print $0,$4/$13}}' $TEMP_OUTPUT_DIR/blast.o > $TEMP_OUTPUT_DIR/good_blast.o`;
#    if ($?) { die "ERROR executing awk, stoping analysis: error code $?\n"}
}