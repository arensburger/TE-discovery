# Feb 2025 This is the top level script to run the TE discovery pipeline. The inputs are a list of possible
# TE sequences (e.g. from RepeatModeler) and a genome. The scripts necessary for running this should contained
# in the directory $PIPELINE_SCRIPTS 
use strict;
use Getopt::Long;
use File::Temp qw(tempfile);

my $PIPELINE_SCRIPTS = "./TE_discovery_scripts"; # location of all scripts necessary for this pipeline
my $INPUT_TE_SEQUENCES; # fasta formated file with input TE sequences
my $INPUT_GENOME; # fasta formated file with genome that input TE sequences came from
my $ANALYSIS_NAME; # name to give to this analysis, can be any string
my $OUTPUT_DIR="./output_files"; # directory to store output files of current analysis (no slash at the end), these can be destroyed when analysis is finished

### VALIDATE INPUTS read and check the inputs
GetOptions(
	't:s'   => \$INPUT_TE_SEQUENCES,
	'g:s'   => \$INPUT_GENOME,
    'n:s'   => \$ANALYSIS_NAME
);
unless (\$INPUT_TE_SEQUENCES and $INPUT_GENOME and $ANALYSIS_NAME) {
	die "usage perl TEdiscovery.pl <-t fasta file with starting TE sequences REQUIRED> <-g fasta file of genome REQUIRED> <-n name for this analysis>\n";
}

### VALIDATE INPUTS check that the pipeline script folder exists
unless (-d $PIPELINE_SCRIPTS) {
    die "ERROR: Scripts directory below (pointed to in script) cannot be found:\n$PIPELINE_SCRIPTS\n";
}

### VALIDATE INPUTS check that the output folder exists
unless (-d $OUTPUT_DIR) {
    die "ERROR: Folder for output files below (pointed to in script) cannot be found:\n$OUTPUT_DIR\n";
}

### PIPELINE STEP reduce redundancy of TE input sequences
my $TE_non_redundant_filename = $OUTPUT_DIR . "/" . $ANALYSIS_NAME . "-TEnonredundant.fa"; # name of the non redundant file
`cd-hit-est -i $INPUT_TE_SEQUENCES -o $TE_non_redundant_filename -d 0 -aS 0.8 -c 0.8 -G 0 -g 1 -b 500 -T 0 -M 0`;
if ($?) { die "Error executing: cd-hit-est, error code $?\n"}