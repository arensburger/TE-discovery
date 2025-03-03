# Feb 2025 - Top-level script for TE discovery pipeline
# Inputs: List of TE sequences (FASTA) and a genome (FASTA)

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use FindBin::libs;  # this sets up that the directory lib will have all modules necessary to run the program 
use TEdiscovery;    # these are the subcripts necessary to run the pipeline

###      CONSTANTS hardware
my $NUM_THREADS = 8; # number of threads to use

### CONSTANTS filtering which proteins to keep for analysis
###      using Goubert et al. as insparation and a tblastn keep hits that have the following characteristics 
my $GENOME_IDENTITY = 80; # per protein, minimum percent identity between protein and genome
my $COVERAGE_RATIO = 0.5; # per protein, minimum ratio of (blast match length) / (query length)
my $COPY_NUMBER = 2; # minimum number of copies that hit different parts of the genome 
my $MIN_DISTANCE = 10000; # if two elements are on the same chromosome, how far they have to be, to be considered different elements

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

### PIPELINE STEP 1 identify proteins that match the genome with high percent identity
my %protein_ids; # holds the id the input proteins that passed the filtering tests as key and the number of copies that passed the test as values
my $step_number = 1;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## Excute the tblastn search
    `tblastn -query $INPUT_PROTEIN_SEQUENCES -db $INPUT_GENOME -outfmt "6 qseqid sseqid sstart send pident length qlen" -out $TEMP_OUTPUT_DIR/tblastn.o -num_threads $NUM_THREADS`;
    if ($?) { die "ERROR executing tblastn, stoping analysis: error code $?\n"}

    ## Inspired by the Goubert et al. protocol, filter elements that 1) have >= 80% identity to genome, 2) have 50% length of the query, 3) are found at multiple locations
    my %candidate_protein; # hash with protein name as key and string with chromosome and middle location of element on that chromsome

    open (INPUT, "$TEMP_OUTPUT_DIR/tblastn.o") or die "ERROR: Cannot open file $TEMP_OUTPUT_DIR/tblastn.o\n";
    while (my $line = <INPUT>) {
        my $gi=0; # boolean, set to zero until the genome identity test is passed
        my $cr=0; # boolean, set to zero until the coverage ratio test is passed
        my $md=1; # boolean, set to one unless mininum distance test is failed

        my @data = split "\t", $line;
        my $middle_position = ($data[3] + $data[2])/2; # position of this element on this chromosome

        if ($data[4] >= $GENOME_IDENTITY) { # test percent id
            $gi=1;
        }
        if ($data[5]/$data[6] >= $COVERAGE_RATIO) { # test for length of the match
            $cr=1;
        }

        # check if current element is close to a recorded one
        if (exists $candidate_protein{$data[0]}) {
            for my $i ( 0 .. $#{ $candidate_protein{$data[0]} } ) {
                my @d2 = split "\t", $candidate_protein{$data[0]}[$i]; # @d2 holds the locus of a previous blast match
                if ($data[1] eq $d2[0])  { # first check that the current and previous locus are on the same chromosome
                    if ((abs($d2[1] - $middle_position)) <= $MIN_DISTANCE) { # second check if they close to each other
                        $md = 0; 
                    }
                }
            }
        }

        if ($gi and $cr and $md) { # if the current protein and locus combination passed all the tests then record it
            push @{ $candidate_protein{$data[0]}}, "$data[1]\t$middle_position";
            $protein_ids{$data[0]} = $#{$candidate_protein{$data[0]}} + 1; # update the hash %protein_ids with the current number of loci that passed the tests
        }      
    }
    close INPUT;

    # Write the resutls into the output file
    my $output_file_name = "$TEMP_OUTPUT_DIR/results_s1.xls";
    open (OUTPUT, '>', $output_file_name) or die "ERROR: Cannot write to $output_file_name\n";
    my $i=0; # counts the number of output lines, to check if it's zero
    foreach my $prot (keys %protein_ids) {
        if ($protein_ids{$prot} >= $COPY_NUMBER) {
            print OUTPUT "$prot\t$protein_ids{$prot}\n";
            $i++;
        }
    }
    close OUTPUT;

    if ($i) {
        print "Finished STEP $step_number and wrote $i identified candidates into file $output_file_name\n";
    }
    else {
        warn "WARNING: STEP $step_number did not result in any identified candiates, no output produced\n";
    }
}