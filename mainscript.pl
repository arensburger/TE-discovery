# Feb 2025 This is the top level script to run the TE discovery pipeline. The inputs are a list of possible
# TE sequences (e.g. from RepeatModeler) and a genome. The scripts necessary for running this should contained
# in the directory $PIPELINE_SCRIPTS 
use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use FindBin::libs;  # this sets up that the directory lib will have all modules necessary to run the program 
use TEdiscovery;    # these are the subcripts necessary to run the pipeline

my $INPUT_TE_SEQUENCES; # fasta formated file with input TE sequences
my $INPUT_GENOME; # fasta formated file with genome that input TE sequences came from
my $ANALYSIS_NAME; # name to give to this analysis, can be any string
my $TEMP_OUTPUT_DIR; # directory where analsysis output files of the analysis are stored, these are not destined to be permananet

### CHECK INPUTS read and check that the inputs have been provided
GetOptions(
	't:s'   => \$INPUT_TE_SEQUENCES,
	'g:s'   => \$INPUT_GENOME,
    'n:s'   => \$ANALYSIS_NAME
);
unless (\$INPUT_TE_SEQUENCES and $INPUT_GENOME and $ANALYSIS_NAME) {
	die "usage perl TEdiscovery.pl <-t fasta file with starting TE sequences REQUIRED> <-g fasta file of genome REQUIRED> <-n name for this analysis>\n";
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

### PIPELINE STEP reduce redundancy of TE input sequences
my $TE_non_redundant_filename = $TEMP_OUTPUT_DIR . "/" . $ANALYSIS_NAME . "-TEnonredundant.fa"; # name of the non redundant file
`cd-hit-est -i $INPUT_TE_SEQUENCES -o $TE_non_redundant_filename -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 0 -M 0`;
if ($?) { die "Error executing: cd-hit-est, error code $?\n"}