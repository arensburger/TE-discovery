# TE Discovery Pipeline - Main Script
# Author: Peter Arensburger
# Date: March 2025
# Description: This script runs the full TE discovery pipeline using input protein sequences and genome data.

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use List::UtilsBy qw(max_by);
use List::Util qw(max);
use Scalar::Util qw(looks_like_number);
use File::Temp qw(tempdir);
use Term::Prompt;
use Term::ANSIColor;

### INPUTs from command line, top level variables are in uppercase
my $INPUT_PROTEIN_SEQUENCES; # fasta formated file with input protein sequences
my $TBLASTN_FILE; # name of the file containing the out put of the tblastn
my $INPUT_GENOME; # fasta formated file with genome that input proteins
my $ANALYSIS_NAME; # base name for the analysis the analysis folder and element folder will be created based on this
my $ANALYSIS_FOLDER; # name of folder to store analysis files into
my $ELEMENT_FOLDER; # directory where individual folders for each element are stored
my $REJECTED_ELEMENTS_FOLDER = "Rejected_elements"; # name of folder that contains files for elements have been reviewed and rejected
my $STEP; # analysis step to perform
my $START_STEP; # analysis step to start at
my $END_STEP; # analysis step to end at,
my $SHOW_HELP; # call for help 

### CHECK INPUTS Read and check that the inputs have been provided
GetOptions(
	'p:s'   => \$INPUT_PROTEIN_SEQUENCES,
	't:s'   => \$TBLASTN_FILE,
	'g:s'   => \$INPUT_GENOME,
    'n:s'   => \$ANALYSIS_NAME,
    's:s'   => \$STEP,
    'h'     => \$SHOW_HELP,
);

## CHECK INPUTS if help was called (this could probably be improved)
if ($SHOW_HELP) {
    print STDERR "Description and input help can be found at https://github.com/arensburger/TE-discovery\n";
    exit;
}

## CHECK INPUTS Validate required input files
unless ($ANALYSIS_NAME and $STEP) {
    die "usage: perl mainscript.pl <-n analysis name REQUIRED> <-s step number REQUIRED> <step specific parameters, use -h for more help>\n";
}

## CHECK INPUTS (NOT STEP SPECIFIC) AND CREATE DIRECTORIES IF NECESSARY
$ANALYSIS_NAME = fixdirname($ANALYSIS_NAME);
$ANALYSIS_FOLDER = $ANALYSIS_NAME . "-analysis";
$ELEMENT_FOLDER = $ANALYSIS_NAME . "-element";
$ANALYSIS_FOLDER = fixdirname($ANALYSIS_FOLDER);
$ELEMENT_FOLDER = fixdirname($ELEMENT_FOLDER);
my $reject_folder_path = $ANALYSIS_FOLDER . "/" . $REJECTED_ELEMENTS_FOLDER;


if (($STEP == 1) or ($STEP == 12)) {
    if (-d $ANALYSIS_FOLDER) { 
        print "WARNING: folder $ANALYSIS_FOLDER already exists, using this folder rather than creating a new one\n";
    }
    else {
        print STDERR "Creating directory $ANALYSIS_FOLDER for storing analysis files\n";
        `mkdir $ANALYSIS_FOLDER`;
        if ($?) { die "ERROR creating directory: error code $?\n"}
    }

    if (-d $ELEMENT_FOLDER) { 
        print "WARNING: folder $ELEMENT_FOLDER already exists, using this folder rather than creating a new one\n";
    }
    else {
        print STDERR "Creating directory $ELEMENT_FOLDER for storing element files\n";
        `mkdir $ELEMENT_FOLDER`;
        if ($?) { die "ERROR creating directory: error code $?\n"}
    }

    if (-d $reject_folder_path) { 
        print "WARNING: folder $reject_folder_path already exists, using this folder rather than creating a new one\n";
    }
    else {
        print STDERR "Creating directory $reject_folder_path for reject elements\n";
        `mkdir $reject_folder_path`;
        if ($?) { die "ERROR creating directory: error code $?\n"}
    }
}
elsif (($STEP == 2) or ($STEP == 3) or ($STEP == 4)) {
    unless ((-d $ANALYSIS_FOLDER) and (-d $ELEMENT_FOLDER) and (-d $reject_folder_path)) {
        die "ERROR: Expecting folders $ANALYSIS_FOLDER, $ELEMENT_FOLDER, and $reject_folder_path would have been created prior to running step $STEP\n"
    }
}
else {
    die "ERROR: don't recognize the step specified with the -s parameter. Stopping analysis.\n";
}

## Create or open files to store analysis parameters, and to store rejected sequences. 
my $analysis_parameters_file_name = "$ANALYSIS_FOLDER/Analysis_run_and_parameters.txt"; # file to record parameters and how the run went
my $rejection_file_name = "$ANALYSIS_FOLDER/Rejected_sequences.txt"; # file to store rejected sequences, and why
open (ANALYSIS,'>>', $analysis_parameters_file_name) or die "ERROR: cannot open file $analysis_parameters_file_name\n"; # create or append to file
open (REJECT, '>>', $rejection_file_name) or die "ERROR, cannot create output file $rejection_file_name\n"; # create or append to file

my $datestring = localtime();
print ANALYSIS "\nSTARTING ANALYSIS on $datestring\n";



### PIPELINE STEP 1 identify proteins that match the genome with parameters specified above under
###     The output is a list of proteins for further analysis recorded in the file $output_file_name

### CONSTANTS applicable to this step only (also record these in the file)
my $GOOD_BLAST_OUTPUT_FILE_NAME = "$ANALYSIS_FOLDER/good_blast.o"; # blast file filtered for good elements according to Goubert
my $BLAST_OUTPUT_FILE_NAME = "$ANALYSIS_FOLDER/tblastn.o"; # default name and location unless a file is provided
my $GENOME_IDENTITY = 80; # IDENTIFYING PROTEINS, per protein, minimum percent identity between protein and genome
my $COVERAGE_RATIO = 0.5; # IDENTIFYING PROTEINS, per protein, minimum ratio of (blast match length) / (query length)
my $COPY_NUMBER = 3; # IDENTIFYING PROTEINS, minimum number of copies that hit different parts of the genome 
my $MIN_DISTANCE = 10000;   # IDENTIFYING PROTEINS, if two elements are on the same chromosome, how far they have to be, to be considered different elements
                            # NOTE: the minimum distance should bigger than the $BLAST_EXTEND variable, to avoid having the same element recorded twice    
my $NUM_THREADS = `nproc --all`;# determine the number of processors on the current machine
if ($?) { print STDERR "WARNING could not determine the number of cores automatically, defaulting to 8\n"; $NUM_THREADS=8}
chomp $NUM_THREADS; 

if (($STEP == 1) or ($STEP == 12)) { # check if this step should be performed or not  
    print STDERR "Working on STEP 1 ...\n";

    ## update the analysis file with what is going on
    print ANALYSIS "Running STEP 1\n";
    print ANALYSIS "\tGENOME_IDENTITY = $GENOME_IDENTITY\n";
    print ANALYSIS "\tCOVERAGE_RATIO = $COVERAGE_RATIO\n";
    print ANALYSIS "\tCOPY_NUMBER = $COPY_NUMBER\n";
    print ANALYSIS "\tMIN_DISTANCE = $MIN_DISTANCE\n";
    print ANALYSIS "\tNUM_THREADS = $NUM_THREADS\n";

    ## VARIABLES, variable for this step
    my %protein_ids; # holds the id the input proteins that passed the filtering tests as key and the number of copies that passed the test as values
    my %rejected_ids; # id's that did not make the cut

    ## Either excute the tblastn search or load the output of a previous run
    if ($TBLASTN_FILE) { # tblastn file was provided
        # sanity check that the file looks ok
        open (INPUT, $TBLASTN_FILE) or die "Cannot open tblastn file $TBLASTN_FILE\n";
        while (my $line = <INPUT>) {
            unless ($line =~ /^\S+\s\S+\s\d+\s\d+\s\S+\s\d+\s\d+\s$/) {
                die "ERROR: tblastn file $TBLASTN_FILE is not formatted as expected at line\n$line";
            }
        }
        close INPUT;

        $BLAST_OUTPUT_FILE_NAME = $TBLASTN_FILE;
        print ANALYSIS "\ttblastn file was provided in file $TBLASTN_FILE\n";
    }
    else {
        # check that all the necessary files have been provided
        unless ($INPUT_PROTEIN_SEQUENCES and $INPUT_GENOME) {
            die ("ERROR, running STEP 1 requires that either 1) both -p and -g parameters are set or 2) that -t is set\n");
        }

        # update the analysis file
        print ANALYSIS "\tInput file: $INPUT_PROTEIN_SEQUENCES\n";
        print ANALYSIS "\tGenome: $INPUT_GENOME\n";

        # Find duplicate, or near duplicate, sequences in the input protein file
        my $protein_file_no_redudants = File::Temp->new(UNLINK => 1, SUFFIX => '.fa' ); # name of file output of cd-hit
        `cd-hit -i $INPUT_PROTEIN_SEQUENCES -o $protein_file_no_redudants -T 0`;
        if ($?) { die "ERROR executing cd-hit: error code $?\n"};

        # identify sequences that are duplicates and updated the analysis files
        my %cluster_number; # holds the sequence name as key and cluster number as value (excluding top sequence)
        my %cluster_topseq; # holds the cluster number as key and reference element as value
        open (INPUT, "$protein_file_no_redudants.clstr") or die "ERROR: Cannot open cluster file $protein_file_no_redudants.clstr\n";
        my $current_cluster_number;
        while (my $line = <INPUT>) { # record all the relevant information from the .clstr file
            if ($line =~ /^>Cluster\s(\d+)/) {
                $current_cluster_number = $1;
            }
            elsif ($line =~ />(\S+)\.\.\.\s\*/) {
                $cluster_topseq{$current_cluster_number}=$1;
            }
            elsif ($line =~ />(\S+)\.\.\.\sat\s/) {
                $cluster_number{$1}=$current_cluster_number;
            }
            else {
                die "ERROR: unexpected line in cluster file $protein_file_no_redudants.clstr\n$line";
            }
        }
        foreach my $dupseq (keys %cluster_number) { # update the Rejected file
            my $datestring = localtime();
            my $topseq = $cluster_topseq{$cluster_number{$dupseq}};
            print REJECT "$datestring $dupseq overlaps with $topseq and was taken out of the analysis at STEP 1 \n";
        }

        # run tblastn
        `tblastn -query $protein_file_no_redudants -db $INPUT_GENOME -outfmt "6 qseqid sseqid sstart send pident length qlen" -out $BLAST_OUTPUT_FILE_NAME -num_threads $NUM_THREADS >&2`;
        if ($?) { die "ERROR executing tblastn, stopping analysis (hint: was the genome formated with makeblastdb?): error code $?\n"}
        print ANALYSIS "\ttblastn was executed, the output is in file $BLAST_OUTPUT_FILE_NAME\n";
    }

    ## Inspired by the Goubert et al. protocol, filter elements that 1) have >= 80% identity to genome, 2) have 50% length of the query, 3) are found at multiple locations
    my %candidate_protein; # hash with protein name as key and string with chromosome and middle location of element on that chromsome
    open (INPUT, "$BLAST_OUTPUT_FILE_NAME") or die "ERROR: Cannot open file $BLAST_OUTPUT_FILE_NAME\n";
    open (OUTPUT, ">", $GOOD_BLAST_OUTPUT_FILE_NAME) or die "ERROR: Cannot create good blast file $GOOD_BLAST_OUTPUT_FILE_NAME\n";
    while (my $line = <INPUT>) {
        my $gi=0; # boolean, set to zero until the genome identity test is passed
        my $cr=0; # boolean, set to zero until the coverage ratio test is passed
        my $md=1; # boolean, set to one unless minimum distance test is failed

        my @data = split "\t", $line;
        $data[0] =~ s/\s//g;
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

        if ($gi and $cr and $md) { # if the current protein and locus combination passed all the tests then record it in %protein_ids, if not record as rejected
            push @{ $candidate_protein{$data[0]}}, "$data[1]\t$middle_position";
            $protein_ids{$data[0]} = $#{$candidate_protein{$data[0]}} + 1; # update the hash %protein_ids with the current number of loci that passed the tests
            print OUTPUT "$line"; # add the current line to the blast file that will be used down the line
        }   
        else {
            $rejected_ids{$data[0]} = "STEP 1\tError code (genome identity/coverage/minimum distance): $gi $cr $md"
        }   
    }
    close INPUT;
    close OUTPUT;

    # updated the rejected file
    foreach my $r (keys %rejected_ids) {
        unless (exists $protein_ids{$r}) { # an id can have multiple blast lines, some acccepted some rejected
            my $datestring = localtime();
            print REJECT "$datestring\t$r\t$rejected_ids{$r}\n";
        }
    }

    # Filter out elements that have too few copy numbers
    # (making this a separate step so the code is more modular, rather then incoroporating it into the next step)
    # also record any elements that were discarded at this point
     foreach my $prot_name (keys %protein_ids) {
        unless ($protein_ids{$prot_name} >= $COPY_NUMBER) {
            my $datestring = localtime();
            print REJECT "$datestring\t$prot_name\tSTEP 1\tError code (BLAST minimum copy number/observed copy number) $COPY_NUMBER $protein_ids{$prot_name}\n";
            delete $protein_ids{$prot_name};
        }
    }
    close OUTPUT;

    # Create individual directories for each element
    my $i=0; # counts the number of output lines, to check if it's zero
    foreach my $prot_name (keys %protein_ids) {
        unless (-d "$ELEMENT_FOLDER/$prot_name") {
            mkdir( "$ELEMENT_FOLDER/$prot_name" ) or die "Couldn't create $ELEMENT_FOLDER/$prot_name directory, $!";
        }
        $i++;
    }

    if ($i) {
        print STDERR "Finished STEP 1, identified $i candidates for further analysis\n";
    }
    else {
        print STDERR "WARNING: STEP 1 did not result in any identified candiates, no output produced\n";
    }
}

### PIPELINE STEP 2 
### For each element that had enough approved blast hits, look for TSD-TIR <---> TIR-TSD combinations 
### CONSTANTS applicable only for STEP 2
my $BLAST_EXTEND = 2000; # IDENTIFYING PROTEINS, number of bp to extend on each side of the blast hit
my $MAX_SEQUENCE_NUMBER = 100; # ALIGNING SEQUENCES maximum number of sequences to consider, to save time
my $CONSLEVEL=0.60; # MAKING CONSENSUS OF SEQUENCES sequence consensus level for consensus
my $MIN_TIR_SIZE = 10; # IDENTIFYING TIR-TSDS smallest allowable size for the TIR
my $TIR_MISMATCHES = 2; # IDENTIFYING TIR-TSDS maximum number of mismatches allowed between two TIRs
my $MIN_PROPORTION_SEQ_WITH_TIR=0.25; #IDENTIFYING TIRs minimum proportion of total elements for a sequence that must contain proper TIRs to be considered a candidate
my $MAX_TIR_PROPORTION=0.75; #IDENTIFYING TIRs how close to the maximum number of tirs do you have to be to qualify as a top TIR
my $MAX_END_PROPORTION=0.75; #IDENTIFYING TIRs how close to maximum proportion of sequences with identical start and stop of tir sequences you can be to a top number
my $MAX_TSD_PROPORTION=0.5; #IDENTIFYING TIRs how close to maximum number of TSDs to qualify as a top TSD sequence
my %EXAMINE_CODES=("1111" => 1, "1101" => 2,  "1110" => 3); # success codes to examine as key and priority as value
my $HIGH_POSITION_CONSENSUS=0.75; # proportion of conservation at an alignment position to call it highly conserved
my $SEARCH_WINDOW_SIZE=20; # how big a window to search on either side of a potential transition postion
my $MAX_GAP_N_AT_POSITION=0.5; # maximum proportion of gaps or N's at an alignment position, if above the position is ignored in this analysis
my $TIR_SEARCH_RANGE = 6; # how many bp away from the alignement transistion to search for TIRs

if (($STEP == 2) or ($STEP == 12)) { # check if this step should be performed or not  
    print STDERR "Working on STEP 2 ...\n";

    ## update the analysis file with what is being done and paramter values
    print ANALYSIS "Running STEP 2\n";
    print ANALYSIS "\tBLAST_EXTEND = $BLAST_EXTEND\n";
    print ANALYSIS "\tMAX_SEQUENCE_NUMBER = $MAX_SEQUENCE_NUMBER\n";
    print ANALYSIS "\tCONSLEVEL = $CONSLEVEL\n";
    print ANALYSIS "\tMIN_TIR_SIZE = $MIN_TIR_SIZE\n";
    print ANALYSIS "\tTIR_MISMATCHES = $TIR_MISMATCHES\n";
    print ANALYSIS "\tMIN_PROPORTION_SEQ_WITH_TIR = $MIN_PROPORTION_SEQ_WITH_TIR\n";
    print ANALYSIS "\tMAX_TIR_PROPORTION = $MAX_TIR_PROPORTION\n";
    print ANALYSIS "\tMAX_END_PROPORTION = $MAX_END_PROPORTION\n";
    print ANALYSIS "\t%EXAMINE_CODES=(\"1111\" => 1, \"1101\" => 2, \"1110\" => 3)\n";
    print ANALYSIS "\tMAX_TSD_PROPORTION = $MAX_TSD_PROPORTION\n";
    print ANALYSIS "\tSEARCH_WINDOW_SIZE = $SEARCH_WINDOW_SIZE\n";
    print ANALYSIS "\tMAX_GAP_N_AT_POSITION = $MAX_GAP_N_AT_POSITION\n";
    print ANALYSIS "\tTIR_SEARCH_RANGE = $TIR_SEARCH_RANGE\n";

    ## check that all the necessary files have been supplied
    # checking that the genome length file is present
    unless ((-f "$INPUT_GENOME.length") and ($INPUT_GENOME)){
        die "ERROR: for this step you need to provide\n1) a fasta formated genome file, using the -g parameter\n2) in the same folder an associated length file generated using the commands below (genome must have been formated using makeblastdb)\n\tsamtools faidx \$genome\n\tawk \'{OFS=\"\\t\"; print \$1,\$2}\' < \$genome.fai > \$genome.length\n";
    }
    print ANALYSIS "\tGenome: $INPUT_GENOME\n";
    print ANALYSIS "\tGenome length file: $INPUT_GENOME.length\n";

    ## VARIABLES, variable for this step
    my @elements; # name of all the elements that will be anlaysed in this step

    ## STEP 2.1
    ## figure out the elements that we're working with 
    ## (for the sake of being modular, redoing this instead of just taking the data from the previous step,
    ## this also ensures that a folder has been created for each element).

    opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
    while (readdir $dh) {
        unless ($_ =~ /^\./) { # prevents reading invisible files or the . and .. filesawk '{OFS="\t"; print $1,$2}' < $genome.fai > $genome.length
            push @elements, $_;
        }
    }
    unless (scalar @elements) { # check that at least one element is present to analyze
        die "ERROR: No elements to analyze were found in the folder $ELEMENT_FOLDER\n";
    }

    ## STEP 2.2
    ## identify TSD-TIR combinations for each element
    my $i=1; # counter of which element is currently being analyzed 
    foreach my $element_name (@elements) {

        # report progess to screen
        my $size = scalar @elements;
        print STDERR "\tProcessing $element_name, number $i of $size\n";
        $i++;

        # create or open the README file for this element
        open (README, ">$ELEMENT_FOLDER/$element_name/$element_name-README.txt") or die "ERROR: Could not open or create README file $ELEMENT_FOLDER/$element_name/$element_name-README.txt\n";

        # STEP 2.2.1
        # for each sequence of this element extend it by $BLAST_EXTEND bps on both sides of the sequence
        my @blastlines = (); # holds all the relevant blast lines for this element
        my $j; # counter of the number of blast lines for this element
        open (INPUT, $GOOD_BLAST_OUTPUT_FILE_NAME) or die "ERROR: cannot open file good blast file $GOOD_BLAST_OUTPUT_FILE_NAME\n";
        while (my $line = <INPUT>) {
            my @data = split ' ', $line;
            if ($data[0] eq $element_name) {
                if ($j < $MAX_SEQUENCE_NUMBER) { # this will ensure that not too many blast results will be used in the analysis, which would slow down the analysis
                    push @blastlines, $line;
                }
                $j++;
            }
        }
        close INPUT;
        if ($j>= $MAX_SEQUENCE_NUMBER) {
            my $datestring = localtime();
            print README "$datestring, The number of BLAST hits ($j) exceeded the maximum for analysis ($MAX_SEQUENCE_NUMBER) analyzing only the first $MAX_SEQUENCE_NUMBER sequences\n";    
        }

        # create the bed file
        open (OUTPUT, '>', "$ELEMENT_FOLDER/$element_name/$element_name.bed") or die "$!\n"; # save the bed file of the original elements that started the analysis
        my $datestring = localtime();
        print README "$datestring, File $element_name.bed contains the positions of all the sequences for element from the BLAST file, this is prior to extension\n";    
 
        foreach my $line (@blastlines) {
            my @data = split ' ', $line;
            if ($data[2] < $data[3]) {
                print OUTPUT "$data[1]\t$data[2]\t$data[3]\t$data[0]\t.\t+\n";
            }
            elsif ($data[2] > $data[3]) {
                print OUTPUT "$data[1]\t$data[3]\t$data[2]\t$data[0]\t.\t-\n";
            }
            else {
                die "ERROR, boundaries of blast cannot be interpreted for line\n$line"
            }
        }
        close OUTPUT;

        # create the files with extended boundaries using bedtools
        my $slopfile = File::Temp->new(UNLINK => 1, SUFFIX => '.slop' );
        `bedtools slop -s -i "$ELEMENT_FOLDER/$element_name/$element_name.bed" -g "$INPUT_GENOME.length" -b $BLAST_EXTEND > $slopfile`;
        if ($?) { die "ERROR executing bedtools: error code $?\n"}

        # for very small chromosomes extending the blast hits can result in duplicate lines in the slop file, removing any duplicate lines
        my $slopfile2 = File::Temp->new(UNLINK => 1, SUFFIX => '.slop' );
        `sort $slopfile | uniq > $slopfile2`;
        if ($?) { die "ERROR removing duplicate lines from slop file: error code $?\n"}

         my $extended_fasta_name = File::Temp->new(UNLINK => 1, SUFFIX => '.fa' ); # name fo the file with extended fasta sequences
        `bedtools getfasta -fi $INPUT_GENOME -fo $extended_fasta_name -bed $slopfile2 -s`;
        if ($?) { die "ERROR executing bedtools: error code $?\n"}

        # 2.2.2 Align the sequences
        my $aligned_sequences_file_name = "$ELEMENT_FOLDER/$element_name/$element_name.maf";
        `mafft --quiet --thread -1 $extended_fasta_name > $aligned_sequences_file_name`;
        if ($?) { die "Error executing mafft, error code $?\n"}
        my $datestring = localtime();
        print README "$datestring, Aligned extended BLAST sequences are in file $element_name.maf\n";
        
        # 2.2.3 determine the highest percentage of agreement on a single nucleotide at each position
        # go through the alignment to record positions that have high agreement on a single nucleotide.

        my %aliseq = fastatohash($aligned_sequences_file_name); # aligned sequences
        my $alignment_length = length($aliseq{(keys %aliseq)[rand keys %aliseq]}); # pick a random sequence to get the length of the alignment (assuming all are the same length)
        my @agreement_proportion; # holds the highest percentage of sequences that agree on one nucleotide at a postion
        my %location_conversion; # holds the position of the alignment without gaps as key and corresponding alignment with gaps as value
        my $consensus_sequence; # this will hold a consensus sequence for the whole alignment, in this consensus "n" means there's too many gaps or n sequences, "N" means there are enough nucleotides but they don't agree

        for (my $i=0; $i<$alignment_length; $i++) { # go through each position in the alignment
            my %chars; # holds the characters found at the current position as key, and abundance as value
            foreach my $taxon (keys %aliseq) {
                my $character = lc(substr($aliseq{$taxon}, $i, 1));
                $chars{$character} += 1;
            }

            my $number_of_sequences_in_alignment = keys %aliseq;            
            # decide if this position has too many N's or gaps
            if ((($chars{"-"}/$number_of_sequences_in_alignment) > $MAX_GAP_N_AT_POSITION) or (($chars{"n"}/$number_of_sequences_in_alignment) > $MAX_GAP_N_AT_POSITION)){
                $consensus_sequence .= "n";
            }
            else { # go here if this position is ok to process further
                my $most_abundant_character = max_by { $chars{$_} } keys %chars;
                if (($most_abundant_character eq "-") or ($most_abundant_character eq "n")) {
                    push @agreement_proportion, 0;
                }
                else {

                    # get the proportion of sequences that agree on the most abundant character, exculding gaps
                    push @agreement_proportion, ($chars{$most_abundant_character}/((keys %aliseq) - $chars{"-"}));                   
                }
                $location_conversion{scalar @agreement_proportion} = $i+1; # update the convertion hash so the correct position can be recorded later

                #update the consensus sequence
                if ($agreement_proportion[-1] >= $CONSLEVEL) {
                    $consensus_sequence .= $most_abundant_character;
                }
                else {
                    $consensus_sequence .= "N"; 
                }
            }
            
        }
        # add the consensus sequence to the alignment
        my $temp_consensus_file = File::Temp->new(UNLINK => 1); # hold the consensus sequence temporarily in this file
        open (OUTPUT, ">", $temp_consensus_file) or die $!;
        print OUTPUT ">consensus\n$consensus_sequence\n";
        close OUTPUT;
        `cat $aligned_sequences_file_name >> $temp_consensus_file`;
        `mv $temp_consensus_file $aligned_sequences_file_name`;

        # 2.2.4 Find location with highest likelihood of being the transition. The methodology here is go through each position of the alignment and compare a
        # window of length $SEARCH_WINDOW_SIZE upstream of that position to another window of the same size downstream. For both windows determine if the alignment 
        # agrees on a single sequence for that window. The pair of upstream and downstream windows that 1) differ the highly from each other in alignment agreement measure, and
        # 2) are closest to the edges of the alignment, are identified as the transition positions.

        my $left_highest_transition_position=0; # position of the most likely transition on the left 
        my $left_highest_transition_number; # highest ratio of conserved positions inside / outside the element
        my $right_highest_transition_position=0; # position of the most likely transition on the left 
        my $right_highest_transition_number; # highest ratio of conserved positions inside / outside the element

        for (my $i=0; $i < scalar @agreement_proportion; $i++) { # go through each position that has a consensus nucleotide
            my $cons_left = 0;  # number of positions in the current window above $HIGH_POSITION_CONSENSUS to the left of the current position (not including it)
            my $cons_right = 0; # number of positions in the current window above $HIGH_POSITION_CONSENSUS to the right of the current position (not including it)
            my $cons_current = 0; # boolean, set to 1 if th current position is above $HIGH_POSITION_CONSENSUS

            for (my $j=$i-$SEARCH_WINDOW_SIZE; $j<$i+$SEARCH_WINDOW_SIZE+1; $j++) { # check positions up and down from current position
                unless ((($i-$SEARCH_WINDOW_SIZE) < 0) or (($i+$SEARCH_WINDOW_SIZE+1) > scalar @agreement_proportion)) { # this will exclude searches outside the bounds of the @agreement_proportion array (i.e. below zero or above the size of the array)
                    if ($agreement_proportion[$j] > $HIGH_POSITION_CONSENSUS) { # true if this is a high consensus position in the window
                        if ($j < $i) {
                            $cons_left++;
                        }
                        elsif ($j > $i) {
                            $cons_right++;
                        }
                        else {
                            $cons_current=1;
                        }
                    }
                }    
            }

            # Update to see if better transition has been found
            if ($cons_current) { # only consider transitions at positions with high agreement
                if ((($cons_right+1) - $cons_left) > $left_highest_transition_number) {
                    $left_highest_transition_position = $location_conversion{$i};
                    $left_highest_transition_number = (($cons_right+1) - $cons_left); # the + 1 is to account that the current position is inside the element
                }
                if ((($cons_left+1) - $cons_right) > $right_highest_transition_number) {
                    $right_highest_transition_position = $location_conversion{$i};
                    $right_highest_transition_number = (($cons_left+1) - $cons_right) # the + 1 is to account that the current position is inside the element
                }
            }
        }
         
        # STEP 2.2.5 
        # identify the TIR and TSD locations around the edges of transitions (if transtions were found)
        if (($left_highest_transition_position and $right_highest_transition_position) and ($left_highest_transition_position < $right_highest_transition_position)) { # only continue if transitions were found positioned correctly
            # go up and down the length of $TIR_SEARCH_RANGE from the transition locations and identify the combinations of TIRs and TSDs. The TIRs are identified 
            # using the sequences with large gaps removed (i.e. %seqrmg), while the TSD are identified on the original sequences (i.e. %seqali)
            my $max_TIR_number; # highest number of TIRs observed for one pair of start and end positions
            my $max_proportion_first_last_bases; # highest number of locations that start and end with the same bases
            my $max_TSD_number; # highest number of intact TSDs
            my @tsd_tir_combinations; # array of possible locations along with various information
            for (my $i=-$TIR_SEARCH_RANGE; $i<=$TIR_SEARCH_RANGE; $i++) {
	            for (my $j=-$TIR_SEARCH_RANGE; $j<=$TIR_SEARCH_RANGE; $j++) {
		            my $number_of_tirs_found=0; # nubmer of sequences that match the TIR criteria
                    my %tir_first_and_last_bases; # first and last set of bases of tir as key and abundance as value
                    my %tsds_found; # keys is TSD type "TA", "2", ... "10" and key is number of TSDs found
                    my %tsds_found2; # keys is TSD type "TA", "2", ... "10" and key is number of TSDs found
                    foreach my $sequence_name (keys %aliseq) {
                        my ($tir1_sequence, $tir2_sequence) = gettir($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j,$MIN_TIR_SIZE, $TIR_MISMATCHES);

                        if ($tir1_sequence) {
                            $number_of_tirs_found += 1;
                            my $bases = substr($tir1_sequence, 0, 3) . substr($tir2_sequence, -3, 3); # recording the first and last 3 bases 
                            $tir_first_and_last_bases{$bases} += 1;
                        }

                        $tsds_found{"TA"} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "TA");
                        $tsds_found{2} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "2");
                        $tsds_found{3} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "3");
                        $tsds_found{4} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "4");
                        $tsds_found{5} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "5");
                        $tsds_found{6} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "6");
                        $tsds_found{7} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "7");
                        $tsds_found{8} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "8");
                        $tsds_found{9} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "9");
                        $tsds_found{10} += gettsd($aliseq{$sequence_name}, $left_highest_transition_position+$i, $right_highest_transition_position+$j, "10");
                    }                    

                    # for this combination of positions, what is the highest proportion of sequences that have a particular TIR (determined only by the first 3 bps.) 
                    my $most_abundant_tir_proportion=0;
                    foreach my $name (sort { $tir_first_and_last_bases{$a} <=> $tir_first_and_last_bases{$b} } keys %tir_first_and_last_bases) {
                        $most_abundant_tir_proportion = $tir_first_and_last_bases{$name}/$number_of_tirs_found;
                    }

                    # record all the data for this combination of $i and $j and update maximum values accross different $i and $j's
                    $max_TIR_number = max($max_TIR_number, $number_of_tirs_found); # assign highest of all the TIR proportions 
                    $max_proportion_first_last_bases = max ($max_proportion_first_last_bases, $most_abundant_tir_proportion);
                    $max_TSD_number = max ($max_TSD_number, $tsds_found{"TA"}, $tsds_found{2}, $tsds_found{3}, $tsds_found{4}, $tsds_found{5}, $tsds_found{6}, $tsds_found{7}, $tsds_found{8}, $tsds_found{9}, $tsds_found{10});
                    my $lp = $left_highest_transition_position+$i;
                    my $rp = $right_highest_transition_position+$j;
                    push @tsd_tir_combinations, "$lp\t$rp\t$number_of_tirs_found\t$most_abundant_tir_proportion\t$tsds_found{\"TA\"}\t$tsds_found{2}\t$tsds_found{3}\t$tsds_found{4}\t$tsds_found{5}\t$tsds_found{6}\t$tsds_found{7}\t$tsds_found{8}\t$tsds_found{9}\t$tsds_found{10}";
                }
            }
            # go through the element and identify those candidate locations that pass the tests for TIR-TSD combinations
            my @successful_candidates; # locations and tsd numbers of candidate locations that the analysis will continue with
            my %failed_candidates; # success codes of the failed candidate and number of candidates as value, used to report failure to the user 
            foreach my $candidate (@tsd_tir_combinations) {
                my $min_prop_seq_wtir=0; # boolean set to zero until passes test for $MIN_PROPORTION_SEQ_WITH_TIR;
                my $top_tir_number=0; # boolean set to zero until passes test for being one of the top tir numbers for these sequences
                my $top_end_proportions=0; # boolean set to zero until passes test for having one of the top proportion of TIR that start and stop with the same sequence
                my $top_number_tsds=0; # boolean set to zero until passes test for having a top number of intact tsds     
                
                my @d = split(" ", $candidate);
                my $number_of_tirs = $d[2];
                my $proportion_of_same_start_stop = $d[3];
                my $number_tsds = max($d[4],$d[5],$d[6],$d[7],$d[8],$d[9],$d[10],$d[11],$d[12],$d[13]);
                # test 1 are the number of sequences with TIR high enough?
                if (($number_of_tirs/(keys %aliseq)) >= $MIN_PROPORTION_SEQ_WITH_TIR) {
                    $min_prop_seq_wtir = 1;
                }

                # test 2 is this candidate have one of the most abundant TIR numbers?
                if (($number_of_tirs) >= ($MAX_TIR_PROPORTION * $max_TIR_number)) {
                    $top_tir_number = 1;
                }

                # test 3 does this candidate have one of the highest proportion of TIRs with identical start and end sequences?
                if (($proportion_of_same_start_stop) >= ($MAX_END_PROPORTION * $max_proportion_first_last_bases)) {
                    $top_end_proportions = 1;
                }

                # test 4 does this candidate have a high number of TSDs?
                if ($number_tsds >= ($MAX_TSD_PROPORTION * $max_TSD_number)) {
                    $top_number_tsds = 1;
                }

                # record the results for this candidate location if the code is one of the acceptables ones for further research
                my $current_success_code = $min_prop_seq_wtir . $top_tir_number . $top_end_proportions . $top_number_tsds;
                if (exists $EXAMINE_CODES{$current_success_code}) {
                   push @successful_candidates, "$current_success_code\t$d[0]\t$d[1]\t$d[2]\t$d[4]\t$d[5]\t$d[6]\t$d[7]\t$d[8]\t$d[9]\t$d[10]\t$d[11]\t$d[12]\t$d[13]";
                }
                else  {
                   $failed_candidates{$current_success_code} += 1; 
                }
            }

            # create the .tirtsd file if successful candidates were found, if not update the user and remove element
            if (scalar @successful_candidates) { # candidates were found
                open (OUTPUT,'>', "$ELEMENT_FOLDER/$element_name/$element_name.tirtsd") or die "Error: cannot create file $ELEMENT_FOLDER/$element_name/$element_name.tirtsd, $!\n";
                print OUTPUT "# success_code\tloc1\tloc2\tnumber_of_tirs\tnumber_tsds_TA_through_10\n";
                foreach my $candidate_line (@successful_candidates) {
                    print OUTPUT "$candidate_line\n";
                }          
                close OUTPUT;

                # update README
                my $datestring = localtime();
                print README "$datestring, File $element_name.tirtsd contains the location information of elements that passed the TIR and TSD codes. The number of failed codes is ";
                foreach my $code (keys %failed_candidates) {
                    print README "$code $failed_candidates{$code}| ";
                }
                print README "\n";
            }
            else { # candidates were not found
                my $datestring = localtime(); 
                print README "$datestring, No locations were found that passed the TIR and TSD codes, stopping analysis here.  The number of failed codes is ";
                foreach my $code (keys %failed_candidates) {
                    print README "$code $failed_candidates{$code}| ";
                }
                print README "\n";

                print REJECT "$datestring\t$element_name\tSTEP 2\tNo lines with acceptable TIR and TSD codes were found\n";
                `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
                if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FOLDER: error code $?\n"}
            }
        } 
        else {
            my $datestring = localtime();
            if ($left_highest_transition_position) {
                print README "$datestring, A left element edge was identified at alignment position $left_highest_transition_position but none on the right, stopping analysis here\n";
                print REJECT "$datestring\t$element_name\tSTEP 2\tA left element edge was identified at alignment position $left_highest_transition_position but none on the right\n";
            }
            elsif ($right_highest_transition_position) {
                print README "$datestring, A right element edge was identified at alignment position $right_highest_transition_position but none on the left, stopping analysis here\n";
                print REJECT "$datestring\t$element_name\tSTEP 2\tA right element edge was identified at alignment position $right_highest_transition_position but none on the left\n";
            }
            else {
                print README "$datestring, No element edge was identified on either left or right, stopping analysis here\n";
                print REJECT "$datestring\t$element_name\tSTEP 2\tNo element edge was identified on either left or right\n";
            }
            `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
            if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FOLDER: error code $?\n"}
        }    
        close README; 
    }    
}

### PIPELINE STEP 3 
### Present the elements to the user for manual review
### CONSTANTS applicable only for STEP 3
my $TIR_bp = 30; # how many bp to display on the TIR side
my $FURTHER_REVIEW_FOLDER_NAME = $ANALYSIS_NAME . "-further_review"; # if the user choose to put elements for further review they will go here

if ($STEP == 3) { # check if this step should be performed or not  
    print STDERR "Working on STEP 3 ...\n";

     ## update the analysis file with what is going on
    print ANALYSIS "Running STEP 3\n";

    my $pkey; # pressed key, used for input from user
    my $elements_left_to_review=1; # number of elements left to review, set to a non-zero value initially so that an initial evaluation will be done
    while ($elements_left_to_review) {

        # determine how many elements are left to review, record a name of an element to review, and
        # record details about what's already been reviewed
        $elements_left_to_review = 0; # recheck that there are still elements to review
        opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
        my %seen_tirs; # hash of arrays that has the unique tsd-tir1-tir2 string as key, and array of with data on this combination
        my $element_name; # if there is an element to review, holds the name of the next element to work on
        my $filename_tirtsd; # name of the current .tirtsd file
        my $filename_maf; # name fo the current .maf 
        
        while (readdir $dh) {
            unless ($_ =~ /^\./) { # prevents reading invisible files or the . and .. files     
                # check if a manual review is already present in the README file for this element
                my $current_element_folder = $ELEMENT_FOLDER . "/" . $_ ; # folder with specific element of interest   
                my $grep_res = `grep "Review 1 result" $current_element_folder/*README.txt`;
                if ($grep_res) { # yes, this element has already been reviewed
                    # recording any TIRs
                    if ($grep_res =~ /TSD\s(\S+),\sTIRs\s(\S+)\sand\s(\S+)/) {
                        my $unique_tsdtir_string = $1 . $2 . $3;
                        $seen_tirs{$unique_tsdtir_string}[0]=$1;
                        $seen_tirs{$unique_tsdtir_string}[1]=$2;
                        $seen_tirs{$unique_tsdtir_string}[2]=$3;
                    }
                }
                else { # no, this element has not been reviewed yet
                    unless ($element_name) { # only set a new element to review if it has not been set yet
                        my $passed_checks = 1; # boolean, set to 0 if checks on the expected files are not passed
                        $filename_tirtsd = $ELEMENT_FOLDER . "/" . $_ . "/" . $_ . ".tirtsd"; # name of the current tirtsd file
                        unless (-e $filename_tirtsd) {
                            warn "WARNING: Cannot review element $_ because there's no .tirtsd file found\n";
                            $passed_checks = 0;
                        }
                        $filename_maf = $ELEMENT_FOLDER . "/" . $_ . "/" . $_ . ".maf";
                        unless (-e$filename_maf) {
                            warn "WARNING: Cannot review element $_ because there's no .maf file found\n";
                            $passed_checks = 0;
                        }
                        if ($passed_checks) { # if the expeced file are there, then set the new element to be reviewed
                            $element_name = $_;
                        }
                    }
                    $elements_left_to_review++;
                }
            }
        }

        if ($elements_left_to_review) {
            
            print color('blue');
            print "\nElement $element_name, $elements_left_to_review left to review\n";
            print color('reset');
            print ANALYSIS "\tReview of element $ELEMENT_FOLDER/$element_name\n";
            `pkill java`; # kill a previous aliview window, this could be dangerous in the long run

            # Variables specific to this section
            my $TIR_b1; # left bound of TIR accepted by user
            my $TIR_b2; # right bound of TIR accpted by user
            my $TSD_size; # size of TSD accepted by user
            my $TSD_type; # if empty then it's a number othwise it's TA
            my $TIR_size; # size of TIR 
            my $manual_left_tir; # if the user does more than one round of manual selection, this will save the coordinate from the previous round, otherwise it's empty
            my $manual_right_tir; # if the user does more than one round of manual selection, this will save the coordinate from the previous round, otherwise it's empty
            my $manual_tsd; # if the user does more than one round of manual selection, this will save the coordinate from the previous round, otherwise it's empty
            my $manual_tir_size=0; # if the user does more than one round of manual selection, this will save the coordinate from the previous round, otherwise it's empty
            
            # check the README file for any previous manual review notes and display them
            open (README, ">>$ELEMENT_FOLDER/$element_name/$element_name-README.txt") or die "ERROR: Could not open or create README file $ELEMENT_FOLDER/$element_name/$element_name-README.txt\n";
            my $prior_notes; 
            while (my $line = <README>) {
                if ($line =~ /Manual\sReview\s1\suser\snote/) {
                    $prior_notes .= $line;
                }
            } 
            print color('bold blue');
            if ($prior_notes) {
                print "\nPRIOR REVIEW NOTES FOR ELEMENT $element_name\n$prior_notes";
            }
            print color('reset');

            # Check if the consensus sequence of this element has a TSD-TIR combination that has been seen before.
            # If such a combination has been seen (and only one) then automatically record this as an element
            # without asking the user. In all other cases let the user decide what to do.
            my %nonaligned_to_aligned; # holds the non-aligned position as key and aligned as value, used for reporting positions
            my %alignment = fastatohash("$ELEMENT_FOLDER/$element_name/$element_name.maf");
            my $consensus_sequence = $alignment{"consensus"};
            my $nogap_consensus; # consensus sequence without gaps
            my %known_tsdtirs; # holds the alignment location of known tirs as key, as value holds array [0] tsd, [1] tir length
            for (my $i=0; $i < length $consensus_sequence; $i++) { # create convertion hash
                my $letter = substr $consensus_sequence, $i, 1; # current letter in the alignment
                unless ($letter =~ /N/i) {
                    $nogap_consensus .= $letter;
                    $nonaligned_to_aligned{length $nogap_consensus} = $i+1;
                }
            }
            foreach my $uniquestring (keys %seen_tirs) { # run through all the previously oberved tirs
                my $tir1 = $seen_tirs{$uniquestring}[1];
                my $tir2 = $seen_tirs{$uniquestring}[2];
                if ($nogap_consensus =~ /$tir1.+$tir2/) {
                    my $l1 = $nonaligned_to_aligned{$-[0] + 1}; # left-most match
                    my $l2 = $nonaligned_to_aligned{$+[0]}; # right-most match
                    $known_tsdtirs{"$l1-$l2"}[0] = $seen_tirs{$uniquestring}[0]; # seen TSD
                    $known_tsdtirs{"$l1-$l2"}[1] = length($tir1); # seen TIR length
                    $known_tsdtirs{"$l1-$l2"}[2] = $tir1; # 
                    $known_tsdtirs{"$l1-$l2"}[3] = $tir2; # 
                }
            }

            # record all the relevant information from the .tirtsd file and check if one or more
            # of the proposed locations match a known tisdtir
            my %locs;   # holds base pair locations as key and array of with all the other values
                        # [0-9] number of found TSDs for sizes TA through 10 
            my %observed_locs;  # just like the %locs but holds the information for previously seen TIRs, this will
            my %TA8_locs; # just like above but with only lines that 8bp or TA tsds, the more common ones
            my $TSDmatch=0; # boolean, only goes to 1 at least one TSD-TIR combination matches an observed one

            open (INPUT, $filename_tirtsd) or die "ERROR: Cannot open file $filename_tirtsd, $!";
            <INPUT>; # skip the header
            while (my $line = <INPUT>) { # reading the individual .tirtsd file
                my @d = split " ", $line;
                if (exists $known_tsdtirs{"$d[1]-$d[2]"}) { # if line match a known tsd tir combination
                    for (my $j=0; $j <= 9; $j++) { # record observed TSD numbers
                        $observed_locs{"$d[1]-$d[2]"}[$j] = $d[$j+4];
                    }
                    # check if there's a non-zero number of TSD with size of the observed TSD
                    if (($known_tsdtirs{"$d[1]-$d[2]"}[0] eq "TA") and ($observed_locs{"$d[1]-$d[2]"}[0] > 0)) {
                        $TSDmatch = 1;
                    }
                    elsif (($known_tsdtirs{"$d[1]-$d[2]"}[0] eq "TTAA") and ($observed_locs{"$d[1]-$d[2]"}[3] > 0)) {
                        $TSDmatch = 1;
                    }
                    elsif ($observed_locs{"$d[1]-$d[2]"}[$known_tsdtirs{"$d[1]-$d[2]"}[0]-1] > 0) {
                        $TSDmatch = 1;
                    }
                }
                elsif (($d[4] > 0) or ($d[11] > 0)) { # this line has a TA or 8 bps TSD
                    for (my $j=0; $j <= 9; $j++) { # record observed TSD numbers
                        $TA8_locs{"$d[1]-$d[2]"}[$j] = $d[$j+4];
                    }
                }
                else { # all other lines
                    for (my $j=0; $j <= 9; $j++) { # record observed TSD numbers
                        $locs{"$d[1]-$d[2]"}[$j] = $d[$j+4];
                    }
                }
            }

            # Decide if the element can be processed automatically or needs to be reviewed by the user
            if (((keys %known_tsdtirs) == 1) and ($TSDmatch)) { # this element can be recorded automatically
                print "\tThis element has a TSD and TIR combination that has been recorded previously, recording and moving to the next element\n";
                foreach my $l (keys %known_tsdtirs) { # there's only one line, but this allows me not to have to look up the name of the key
                    my $tsd = $known_tsdtirs{$l}[0];
                    my $tir1 = $known_tsdtirs{$l}[2];
                    my $tir2 = $known_tsdtirs{$l}[3];
                    print README "$datestring, Automated Review 1 result: This is an element, TSD $tsd, TIRs $tir1 and $tir2\n";
                }
            }
            else { # this element need manual review
                # setup and display menu 1
                my %alignment_sequences = fastatohash($filename_maf); # load the existing alignment
                my @menu1_items; # hold the text of menu 1 choices
                push @menu1_items, "Save current progress and quit the program"; # item 0
                push @menu1_items, "Done reviewing this element and move on the next one"; # item 1
                push @menu1_items, "View the whole sequence alignment"; # item 2
                push @menu1_items, "Make a note in the README file"; # item 3
                push @menu1_items, "Update the README to say this is not an element and quit this element"; # item 4
                push @menu1_items, "Set this element aside for later review"; # item 5
                push @menu1_items, "Manually enter TIRs and TSD"; # item 6
                foreach my $loc (keys %observed_locs) {
                    my @pos = split "-", $loc; # get the TIR end positions
                    my ($TIR_number, $average_TIR_length) = average_tir_number_and_length($pos[0], $pos[1], $MIN_TIR_SIZE, $TIR_MISMATCHES, %alignment_sequences);
                    push @menu1_items, "$loc, $average_TIR_length, $observed_locs{$loc}[0]-$observed_locs{$loc}[1]-$observed_locs{$loc}[2]-$observed_locs{$loc}[3]-$observed_locs{$loc}[4]-$observed_locs{$loc}[5]-$observed_locs{$loc}[6]-$observed_locs{$loc}[7]-$observed_locs{$loc}[8]-$observed_locs{$loc}[9] Seen TIRs with $known_tsdtirs{$loc}[0] TSDs"; 
                }
                foreach my $loc (keys %TA8_locs) {
                    my @pos = split "-", $loc; # get the TIR end positions
                    my ($TIR_number, $average_TIR_length) = average_tir_number_and_length($pos[0], $pos[1], $MIN_TIR_SIZE, $TIR_MISMATCHES, %alignment_sequences);
                    push @menu1_items, "$loc, $average_TIR_length, $TA8_locs{$loc}[0]-$TA8_locs{$loc}[1]-$TA8_locs{$loc}[2]-$TA8_locs{$loc}[3]-$TA8_locs{$loc}[4]-$TA8_locs{$loc}[5]-$TA8_locs{$loc}[6]-$TA8_locs{$loc}[7]-$TA8_locs{$loc}[8]-$TA8_locs{$loc}[9]"; 
                }
                foreach my $loc (keys %locs) {
                    my @pos = split "-", $loc; # get the TIR end positions
                    my ($TIR_number, $average_TIR_length) = average_tir_number_and_length($pos[0], $pos[1], $MIN_TIR_SIZE, $TIR_MISMATCHES, %alignment_sequences);
                    push @menu1_items, "$loc, $average_TIR_length, $locs{$loc}[0]-$locs{$loc}[1]-$locs{$loc}[2]-$locs{$loc}[3]-$locs{$loc}[4]-$locs{$loc}[5]-$locs{$loc}[6]-$locs{$loc}[7]-$locs{$loc}[8]-$locs{$loc}[9]"; 
                }

                my $menu1 = 1; # boolean, set to one until the user is done with menu 1
                my $move_to_menu2 = 0; # boolean, set to zero until the user is ready to move on to menu 2
                while ($menu1) {
                    my $menu1_choice = prompt('m', {
                        title => "MENU1 of $element_name\n",
                        prompt => 'What would you like to do?',
                        display_base => 0,
                        return_base => 0,
                        accept_multiple_selections => 0,
                    items  => [@menu1_items],
                    separator => '[,/\s]',
                    },'', 2);
        
                    # process answer to menu 1
                    if ($menu1_choice == 0) {
                        `pkill java`; # kill a previous aliview window, this could be dangerous in the long run
                        close README;
                        close ANALYSIS;
                        close REJECT;
                        exit;
                    }
                    if ($menu1_choice == 1) { # the user is done with this element
                        $menu1 = 0;
                    }
                    elsif ($menu1_choice == 2) { # the user wants to see the original alignment
                        `aliview $ELEMENT_FOLDER/$element_name/$element_name.maf`;
                        if ($?) { die "ERROR: Could not open program aliview: error code $?\n"}
                    }
                    elsif ($menu1_choice == 3) { # the user wants to update the README file or set this element aside
                        
                        # update the README file with a new line and open the editor
                        my $datestring = localtime(); 
                        `echo $datestring, Manual Review 1 user note: >> $ELEMENT_FOLDER/$element_name/$element_name-README.txt`;
                        if ($?) { die "ERROR: could not add line to README file $ELEMENT_FOLDER/$element_name/$element_name-README.txt: error code $?\n"}
                        `gnome-text-editor $ELEMENT_FOLDER/$element_name/$element_name-README.txt`;
                        if ($?) { die "ERROR: could not open for editing the file $ELEMENT_FOLDER/$element_name/$element_name-README.txt: error code $?\n"}
                    }
                    elsif ($menu1_choice == 4) { # the user wants to label this an not an element
                        my $datestring = localtime(); 
                        print README "$datestring, Manual Review 1 result: This is not an element\n";
                        print REJECT "$datestring\t$element_name\tSTEP 3\tManual review of TSD and TIRs determined this is not an element\n";
                        `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
                        $menu1 = 0;
                    }
                    elsif ($menu1_choice == 5) { # the user wants to set this element aside
                        unless (-d $FURTHER_REVIEW_FOLDER_NAME) { # create the folder if necessary
                            `mkdir $FURTHER_REVIEW_FOLDER_NAME`;
                            if ($?) { die "ERROR: could make folder $FURTHER_REVIEW_FOLDER_NAME: error code $?\n"}
                        }
                        print README "$datestring, in STEP 3 manual review 1, reviewer moved this element to the folder $FURTHER_REVIEW_FOLDER_NAME for further review\n";
                        `mv $ELEMENT_FOLDER/$element_name $FURTHER_REVIEW_FOLDER_NAME`;
                        if ($?) { die "ERROR: could not move folder $ELEMENT_FOLDER/$element_name to folder $FURTHER_REVIEW_FOLDER_NAME: error code $?\n"}
                        print "\tMoved the element to the folder $FURTHER_REVIEW_FOLDER_NAME\n";
                        $menu1 = 0;
                    }
                    elsif ($menu1_choice == 6) { # the user wants to manually enter the coordinates
                        $TIR_b1 = prompt('n', "Left alignement coordinate:", '', $manual_left_tir );
                        $TIR_b2 = prompt('n', "Right alignement coordinate:", '', $manual_right_tir);
                        $TSD_size = prompt('n', "TSD size (enter 1 for \"TA\" or 0 for no TSDs):", '', $manual_tsd);
                        $TIR_size = prompt('n', "TIR size (enter 0 for the script to find the most likely size):", '', $manual_tir_size);

                        # set this manual selections for next round
                        $manual_left_tir = $TIR_b1;
                        $manual_right_tir = $TIR_b2;
                        $manual_tsd = $TSD_size;
                        $manual_tir_size = $TIR_size;
                        
                        if ($TSD_size == 1) {
                            $TSD_size = 2;
                            $TSD_type = "TA";
                        }
                        elsif ($TSD_size == 0) {
                            $TSD_size = 0;
                            $TSD_type = "NA";
                        }
                        else {
                            $TSD_type = $TSD_size;
                        }
                        $move_to_menu2 = 1;
                    }
                    else { # This means that the user selected a preset number
                        $move_to_menu2 = 1;
                        if ($menu1_items[$menu1_choice] =~ /^(\d+)-(\d+),\s(\d+),\s(\S+)/) {
                            $TIR_b1 = $1; 
                            $TIR_b2 = $2; 
                            $TIR_size = $3;
                            
                            # get the sizes of the TSDs and their position in the list to decide the best TSD size to default to
                            my @d = split "-", $4;
                            my @sorted_indices = sort { $d[$b] <=> $d[$a] } 0 .. $#d;
                            my $default_tsd = $sorted_indices[0] + 1;
                            if (($default_tsd == 8) and ($d[0] > (5*$d[7]))) { # special case where TA tsds way outnumber 8bp. TSD, likely a case of the 8bp. tsd being TATATATATA so favoring a TA tsd
                                $default_tsd = 1;
                            }
                            elsif ($d[7] > 0.5*($d[$sorted_indices[0]])) { # special case where the 8bp. tsds are larger than half the largest TSD, then go with 8bp.
                                $default_tsd = 8;
                            }
                            # ask what TSD size to use
                            $TSD_size = prompt('n', "What TSD size should be displayed? (enter 1 for \"TA\" or 0 for no TSDs)", '', $default_tsd);
                            if ($TSD_size == 1) {
                                $TSD_size = 2;
                                $TSD_type = "TA";
                            }
                            elsif ($TSD_size == 0) {
                                $TSD_size = 0;
                                $TSD_type = "NA";
                            }
                            else {
                            $TSD_type = $TSD_size; 
                            }
                            $move_to_menu2 = 1;
                        }
                        else {
                            die "ERROR: Cannot parse selection $menu1_items[$menu1_choice]\n";
                        }
                    }

                    if (($menu1) and ($move_to_menu2)){ # only continue if the user has not elected to quit menu 1 and is ready to move on
                    
                        # The TIRs and TSDs location have now been selected, next create an alignment to display these to the user  

                        # Find and record the conensus sequence, this will be necessary to properly display the TIR sequences
                        my $consensus_sequence;
                        foreach my $seq_name (keys %alignment_sequences) {
                            if ($seq_name =~ /consensus/) { 
                                $consensus_sequence = $alignment_sequences{$seq_name};
                            }   
                        }
                        unless ($consensus_sequence) { die "ERROR: No conensus sequence found, need this to properly display the TIRs\n"}

                        my $temp_alignment_file = "/tmp/$element_name-$TIR_b1-$TIR_b2-$TSD_size.fa"; # tried using perl temporary file system, but aliview will not open those
                        my %TIR_sequences; # element name as key and [0] left TIR sequences displayed [1] right TIR sequences displayed. This will be used in the next menus to display TIR sequences
                        open (OUTPUT, ">", $temp_alignment_file) or die "Cannot create temporary alignment file $temp_alignment_file\n";

                        foreach my $seq_name (keys %alignment_sequences) {
                            if ($seq_name =~ /consensus/) { 
                                $consensus_sequence = $alignment_sequences{$seq_name};
                            }
                            else {# avoid the line with the consensus sequence
                                ## left side sequences
                                my $left_whole_seq = substr($alignment_sequences{$seq_name}, 0, $TIR_b1-1);
                                $left_whole_seq =~ s/-//g; #remove gaps
                                # if there are no or few sequences, replace left TSD with space symbols
                                if ((length $left_whole_seq) < $TSD_size) {
                                    $left_whole_seq = "";
                                    for (my $i=0; $i<=$TSD_size; $i++) {
                                        $left_whole_seq .= "s";
                                    }
                                }
                                my $left_tsd = substr($left_whole_seq, -$TSD_size, $TSD_size);

                                # get the sequence of the TIR, ignoring positions with no consensus
                                my $i=0;
                                my $left_tir_seq;
                                while ((length $left_tir_seq) < $TIR_bp) {
                                    unless ((substr $consensus_sequence, $TIR_b1-1+$i, 1) =~ /n/) {
                                        $left_tir_seq .= substr($alignment_sequences{$seq_name}, $TIR_b1-1+$i, 1);
                                    }
                                    $i++;
                                    if ($i > length $alignment_sequences{$seq_name}) { die "ERROR: Cannot display TIR for sequence $seq_name\n"} # a reality check in case things go south
                                }

                                ## right side sequences
                                my $right_whole_seq = substr($alignment_sequences{$seq_name}, $TIR_b2);
                                $right_whole_seq =~ s/-//g; #remove gaps
                                # if there are no or few sequences, replace right TSD with space symbols
                                if ((length $right_whole_seq) < $TSD_size) {
                                    $right_whole_seq = "";
                                    for (my $i=0; $i<=$TSD_size; $i++) {
                                        $right_whole_seq .= "s";
                                    }
                                }
                                my $right_tsd = substr($right_whole_seq, 0, $TSD_size);
                                my $i=0;
                                my $right_tir_seq;
                                while ((length $right_tir_seq) < $TIR_bp) {
                                    unless ((substr $consensus_sequence, $TIR_b2-$i-1, 1) =~ /n/i) {
                                        $right_tir_seq .= substr($alignment_sequences{$seq_name}, $TIR_b2-$i-1, 1);
                                    }
                                    $i++;
                                    if ($TIR_b2-$i < 0) { die "ERROR: Cannot display TIR for sequence $seq_name\n"} # a reality check in case things go south
                                }
                                $right_tir_seq = reverse $right_tir_seq; # necessary because sequences were added from the TIR end backward

                                ## print the sequences after checking if there's anything to print
                                my $test_tir1 = $left_tir_seq;
                                my $test_tir2 = $right_tir_seq;
                                $test_tir1 =~ s/-//g; # removing all the gaps to see if there's anything left after removal
                                $test_tir2 =~ s/-//g;
                                if (($test_tir1) or ($test_tir2)) { # only print if there's something in the TIR sequecences
                                    if ($left_tsd eq $right_tsd) { # if the TSDs are the same (and they are not just S's) then add it to the title
                                        unless (($left_tsd =~ /s/) or ($right_tsd =~ /s/)) {
                                            $seq_name .= "-identicalTSDs";
                                        }
                                    }
                                    print OUTPUT ">$seq_name\n", $left_tsd, "sss", $left_tir_seq, "ssssssssssssssssssss", $right_tir_seq, "sss", $right_tsd, "\n";
                                    $TIR_sequences{$seq_name}[0] = $left_tir_seq; # store the sequences for MENU2
                                    $TIR_sequences{$seq_name}[1] = $right_tir_seq; # store the sequences for MENU2
                                }
                            }
                        }
                        close OUTPUT;

                        # MENU 2 display the alignement to the user and ask for evaluation
                        `aliview $temp_alignment_file`;
                        if ($?) { die "Error executing: aliview $temp_alignment_file, error code $?\n"}
                        my $menu2 = 1; # boolean, set to 1 until the user is done with menu 2

                        # if the used asked the script to figure out the TIR size then it will be set to zero
                        # figure out the average size here
                        if ($TIR_size ==0) {
                            my $TIR_number;
                            ($TIR_number, $TIR_size) = average_tir_number_and_length($TIR_b1, $TIR_b2, $MIN_TIR_SIZE, $TIR_MISMATCHES, %alignment_sequences)
                        }

                        # figure out the TIR sequences
                        # don't use the consensus sequence for this, because that contains n characters
                        my $TIR1_sequence;
                        my $TIR2_sequence; 
                        for (my $i=0; $i <$TIR_size; $i++) {
                            my %left_char_abundance; # holds the character as key and abundance as value
                            my %right_char_abundance;
                            foreach my $seqname (keys %TIR_sequences) {
                                my $left_character = substr($TIR_sequences{$seqname}[0], $i, 1); # current left character
                                my $right_character = substr($TIR_sequences{$seqname}[1], length($TIR_sequences{$seqname}[1])-$TIR_size+$i, 1); # current right character
                                unless ($left_character eq "-") {
                                    $left_char_abundance{$left_character} += 1;
                                }
                                unless ($right_character eq "-") {
                                    $right_char_abundance{$right_character} += 1;
                                }
                            }
                            $TIR1_sequence .= max_by { $left_char_abundance{$_} } keys %left_char_abundance;
                            $TIR2_sequence .= max_by { $right_char_abundance{$_} } keys %right_char_abundance;
                        }
                                            
                        # setup and display menu 2, testing to see if this combination of TIRs has been seen before
                        my @menu2_items; # hold the text of menu 2 choices
                        push @menu2_items, "Go back to the previous menu";

                        # testing for this TIR having been seen before, change the option if it has
                        my ($tirmatch, $other_tir1, $other_tir2, $other_tsd) = compare_tirs($TIR1_sequence, $TIR2_sequence, %seen_tirs);

                        # if the TIR_size is N/A no TIR sequences will be reported, this will match the menu to that
                        if ($TIR1_sequence and $TIR2_sequence and ($TSD_type eq "NA")) { 
                            push @menu2_items, "Update the README to record TIRs $TIR1_sequence and $TIR2_sequence";
                        }
                        elsif ($TIR1_sequence and $TIR2_sequence) {
                            push @menu2_items, "Update the README to say this is an element with TSDs of type $TSD_type bp.\n\tand TIRs $TIR1_sequence and $TIR2_sequence";
                        }
                        else {
                            push @menu2_items, "No TIRs were found, update the README to say this is not an element";
                        }
                        
                        push @menu2_items, "Update the README to say this is not an element";
                        while ($menu2) {  #keep displaying until the user ready to leave

                            # if applicable tell the user this TIRs been seen before
                            print color('bold blue');
                            if (abs($tirmatch)==2) {
                                print "\nNOTE: These exact TIRs have been recorded before\n";
                            }
                            if (abs($tirmatch)==1) {
                                print "\nNOTE: Very similar TIRs have been recorded before\n";
                                print "PRIOR TIRs: $other_tir1  $other_tir2\n";
                            }
                            print color('reset');
                            
                            # display menu 2
                            my $menu2_choice = prompt('m', {
                                title => "MENU2 of $element_name\n",
                                prompt => 'What would you like to do?',
                                display_base => 0,
                                return_base => 0,
                                accept_multiple_selections => 0,
                            items  => [@menu2_items],
                            },'', 0);

                            # process the user's choice1
                            if ($menu2_choice == 0) { # user wants to go back to menu 1 
                                $menu2 = 0;
                                $move_to_menu2 = 0;
                            } 
                            elsif ($menu2_choice == 1) { # user wants to report this an element as it is
                                my $datestring = localtime(); 
                                if ($TIR1_sequence and $TIR2_sequence) { # this prevents printing if no TIRs were found
                                    if ($TSD_type eq "TA") {
                                        print README "$datestring, Manual Review 1 result: This is an element, TSD TA, TIRs $TIR1_sequence and $TIR2_sequence\n";
                                    }
                                    elsif ($TSD_type eq "NA") {
                                        print README "$datestring, Manual Review 1 result: TSD NA, TIRs $TIR1_sequence and $TIR2_sequence\n";               
                                    }
                                    else {
                                        print README "$datestring, Manual Review 1 result: This is an element, TSD $TSD_type, TIRs $TIR1_sequence and $TIR2_sequence\n";
                                    }
#                                    $prior_review = 1;
                                }
                                else { # no TIRS were found
                                    print README "$datestring, Manual Review 1 result: This is not an element\n";
                                    print REJECT "$datestring\t$element_name\tSTEP 3\tManual review of TSD and TIRs determined this is not an element\n";
                                    `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
                                    $menu2 = 0;
                                    $menu1 = 0
                                }
                                $menu2 = 0;
                            }
                            elsif ($menu2_choice == 2) {  # user wants to report this as not an element
                                my $datestring = localtime(); 
                                print README "$datestring, Manual Review 1 result: This is not an element\n";
                                print REJECT "$datestring\t$element_name\tSTEP 3\tManual review of TSD and TIRs determined this is not an element\n";
                                `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
                                if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FOLDER: error code $?\n"}
                                $menu2 = 0;
                            }
                        }
                    }            
                }  
            }
            close README;         
        }
        else {
            print "No elements left that need review\n";
        }
    }
}

### PIPELINE STEP 4 
### Identify most likely transposase ORF based on identified TIR sequences

# need to fix to read the new README file names

if ($STEP == 4) { # check if this step should be performed or not  
    print STDERR "Working on STEP 4 ...\n";

    ## Constant for this step
    my $MAX_ELEMENT_SIZE = 5000; # maximum element size
    my $PROTEIN_TBLASTN_CUTOFF = 1e-6; # E-value cutoff for initial protein to match element nucleotide sequence
    my $CONSENSUS_REMOVAL_THRESHOLD = 0.75; # minium number of taxa with a definite nucleotide at a position to make a consensus at this position
    my $CONSENSUS_LEVEL = 0.5; # proportion of agreeing nucleotides to call a consensus 
    
    ## update the analysis file with what is going on
    print ANALYSIS "Running STEP 4\n";

    # making sure all the required information has been provided
    unless ($INPUT_GENOME and $INPUT_PROTEIN_SEQUENCES){
        die "ERROR: for this step you need to provide two pieces of information:
             1) a fasta formated genome file, using the -g parameter
             2) the input protein file using the -p parameter\n";
    }
    print ANALYSIS "\tGenome: $INPUT_GENOME\n";
    print ANALYSIS "\tInput protein file: $INPUT_PROTEIN_SEQUENCES\n";
    print ANALYSIS "\tMAX_ELEMENT_SIZE = $MAX_ELEMENT_SIZE\n";
    print ANALYSIS "\tPROTEIN_TBLASTN_CUTOFF = $PROTEIN_TBLASTN_CUTOFF\n";
    print ANALYSIS "\tCONSENSUS_REMOVAL_THRESHOLD = $CONSENSUS_REMOVAL_THRESHOLD\n";
    print ANALYSIS "\tCONSENSUS_LEVEL = $CONSENSUS_LEVEL\n";

    ## Read the README files and identify TIRs sequences and TSD
    
    # load the element information from the README files
    my %reviewed_tsdtirs;  # holds the element name as key and information on the element ends as array of values. 
                    # specifically [0] = TSD size or type, [1] = TIR1 sequence, [2] = TIR2 sequence
                    # This will record the information in the last line of the README file that has the right format
    opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
    while (readdir $dh) {
        unless ($_ =~ /^\./) {
            if (-e "$ELEMENT_FOLDER/$_/$_-README.txt") {
                open (INPUT, "$ELEMENT_FOLDER/$_/$_-README.txt") or die "ERROR: Cannot open file $!";
                while (my $line = <INPUT>) {
                    if ($line =~ /This\sis\san\selement,\sTSD\s(\S+),\sTIRs\s(\S+)\sand\s(\S+)/) {
                        $reviewed_tsdtirs{$_}[0] = $1; # TSD 
                        $reviewed_tsdtirs{$_}[1] = $2; # first TIR sequence
                        $reviewed_tsdtirs{$_}[2] = $3; # second TIR sequence
                    }
                }
                close INPUT; 
            }
            else {
                print STDERR "WARNING: No README.txt file was found in folder $ELEMENT_FOLDER/$_\n";
            }
            unless (exists $reviewed_tsdtirs{$_}) {
                print STDERR "WARNING: No appropriate line with TIR sequences were found in the README file for element $_\n";
            }
        }
    }
    unless (keys %reviewed_tsdtirs) {
        print STDERR "WARNING: No elements where found for processing in this dataset\n";
    }

    # Cluster the TIR sequences
    my %clustering_info;   # holds a unique cluster name (i.e. cluster1, cluster2, etc.) as key and as value
                    # an array with [0] string with names of all proteins in this cluster
                    # [1] name of the protein with the shortest TIRs [2] length of shortest TIR [3] last TSD observed [4] 0 if all
                    # TSDs are the same 1 if they are not
    my $tircluster_input_file = File::Temp->new(UNLINK => 1); # clustering program input file
    my $tircluster_output_file = File::Temp->new(UNLINK => 1); # clustering program output file
    open (TIRCLUSTER, ">", $tircluster_input_file) or die "Cannot create temporary file $tircluster_input_file\n";
    foreach my $protein_name (keys %reviewed_tsdtirs) {
        my $seq = "$reviewed_tsdtirs{$protein_name}[1]$reviewed_tsdtirs{$protein_name}[2]";
        print TIRCLUSTER ">$protein_name\n$seq\n";
    }
    close TIRCLUSTER;
    `cd-hit-est -i $tircluster_input_file -o $tircluster_output_file`; # run clustering
    if ($?) { die "ERROR executing cd-hit-est: error code $?\n" }

    # interpret the clustering .clstr file and populate the %clustering_info hash 
    my $cluster_name; # current element name;
    open (INPUT, "$tircluster_output_file.clstr") or die "Cannot open cd-hit output file $tircluster_output_file.clstr\n";
    while (my $line = <INPUT>) {
        if ($line =~ />Cluster\s(\d+)/) {
            $cluster_name = "element$1";
        }
        elsif ($line =~ /^\d+\s+(\d+)nt,\s>(\S+)\.\.\./) {
            my $length = $1;
            my $protein_name = $2;

            # add the protein name to the list 
            $clustering_info{$cluster_name}[0] .= $protein_name . " ";

            # record the name of the shortest TIRs in each cluster
            if ((exists $clustering_info{$cluster_name}[2]) and ($length < $clustering_info{$cluster_name}[2])) {
                $clustering_info{$cluster_name}[1] = $protein_name;
                $clustering_info{$cluster_name}[2] = $length;
            }
            unless (exists $clustering_info{$cluster_name}[2]) {
                $clustering_info{$cluster_name}[1] = $protein_name;
                $clustering_info{$cluster_name}[2] = $length;
            }

            # check that the TSDs are the same
            if (exists $clustering_info{$cluster_name}[3]) { # a TSD has been observed previously
                if ($clustering_info{$cluster_name}[3] ne $reviewed_tsdtirs{$protein_name}[0]) {
                    $clustering_info{$cluster_name}[4] = 1;
                }
            }
            else {
                $clustering_info{$cluster_name}[3] = $reviewed_tsdtirs{$protein_name}[0];
                $clustering_info{$cluster_name}[4] = 0;
            }
        }    
        else {
            die "ERROR: Unexpected line in file $tircluster_output_file.clstr\n$line\n";
        }
    }
    # Warn the user of any elements where the TSD sizes do not all agree
    foreach my $key (keys %clustering_info) {
        if ($clustering_info{$key}[4] == 1) {
        warn "WARNING: Clustered elements $clustering_info{$key}[0] do not all have the same TSDs, these elements will be ignored\n";   
        }
    }


    ## Using tTIR sequences for each cluster, identify possible element sequences, then match these to the original protein sequence.
    ## Report potential complete elements.
    my %proteins = fastatohash($INPUT_PROTEIN_SEQUENCES); # this holds the sequence of all the original protein sequences
    my %genome = fastatohash($INPUT_GENOME); # load the genome into memory
    foreach my $cluster_name (keys %clustering_info) { # go through the clusters individually
      
        # get the TIR sequences for this cluster
        my $tir1_seq = lc($reviewed_tsdtirs{$clustering_info{$cluster_name}[1]}[1]);
        my $tir2_seq = lc($reviewed_tsdtirs{$clustering_info{$cluster_name}[1]}[2]);
        my $TIRs_set = 1; # boolean, set to 1 if TIR sequence are ok, or zero if not
        if (($tir1_seq =~/n/i) or ($tir2_seq =~ /n/i)) { # if the tir sequences include "n" characters warn the user and stop processing
            print STDERR "WARNING: The TIR sequences $tir1_seq and $tir2_seq for element(s) $clustering_info{$cluster_name}[0] contain one or more n characters, this is a problem for finding this element in the genome. Ignoring the element(s).\n";
            $TIRs_set = 0;
        }

        # get the TSD length for this cluster
        my $TSD_length = lc($reviewed_tsdtirs{$clustering_info{$cluster_name}[1]}[0]); # all lower case to avoid confusion
        if ($TSD_length eq "ta") {
            $TSD_length = 2;
        }
        elsif ($TSD_length eq "tata") {
            $TSD_length = 4;
        }
        # identify the nucleotide sequences between TIR locations
        if ($TIRs_set) { # only continue if the TIR sequence are ok
            my %element_sequences; # holds the genomic sequence with intact tirs as well as tsd sequences, it's a hash to avoid duplications
            foreach my $chr (keys %genome) { # go through each genome subsection (calling it "chr" here)
                # Identify all the element sequences (including TSDs) in the forward orientation. Converting everything to lower case to avoid confusion with cases.
                # Also providing the name of the chrososome and orientation so that the position of all the elements can recorded.
                my %fw_element_sequences = identify_element_sequence(lc($genome{$chr}), $tir1_seq, $tir2_seq, $MAX_ELEMENT_SIZE, $chr, $TSD_length, "+"); # look for TIRs on the + strand
                %element_sequences = (%element_sequences, %fw_element_sequences); # add elements found on the + strand to %element_sequences
                unless ($tir1_seq eq (rc($tir2_seq))) { # only look on the other strand if the TIRs are not symetrical, symetrical TIR will already have been found
                    my %rc_element_sequences = identify_element_sequence(lc($genome{$chr}), $tir1_seq, $tir2_seq, $MAX_ELEMENT_SIZE, $chr, $TSD_length, "-");
                    %element_sequences = (%element_sequences, %rc_element_sequences); # add any new elements to the hash %element_sequences
                }
            }  

            # Continue the analysis if complete elements have been found
            if (keys %element_sequences) { 
                # make a blast database using the file that contains the nucleotide sequences of the element
                my $database_input_file = File::Temp->new(UNLINK => 1); 
                open (OUTPUT, ">", $database_input_file) or die "$!";
                foreach my $seq (keys %element_sequences) {
                    print OUTPUT ">$seq\n$element_sequences{$seq}\n";
                }   
                close OUTPUT;

                my $tblastn_database_name = File::Temp->new(UNLINK => 1); # file name for the tblastn database name
                `makeblastdb -in $database_input_file -dbtype nucl -out $tblastn_database_name`;
                if ($?) { die "ERROR executing makeblastdb: error code $?\n"}
# modify this part to do a blast with all the proteins involved
print "$clustering_info{$cluster_name}[0]\n"; exit;
                # make a file with the sequence of the original protein used to find the the current element
                my $protein_file = File::Temp->new(UNLINK => 1);
                open (OUTPUT, ">", $protein_file) or die "$!";
                print OUTPUT ">protein\n$proteins{$cluster_name}\n";
                close OUTPUT;

                # run tblastn
                my $tblastn_output = File::Temp->new(UNLINK => 1); 
                `tblastn -query $protein_file -db $tblastn_database_name -out $tblastn_output -outfmt "6 sseqid evalue sstart send"`;
                if ($?) { die "ERROR executing tblastn: error code $?\n"}   
                # interpret the tblastn results, identify element sequences that have low e-value and report them in the same orientation 
                # as the protein 
                my %complete_elements_sequences; # this will hold the identified element sequences as value, and genomic location as value
                open (INPUT, $tblastn_output) or die "$!";           
                while (my $line = <INPUT>) {
                    my @data = split " ", $line;
                    if ($data[1] < $PROTEIN_TBLASTN_CUTOFF) { # only continue if tblastn E value is low enough
                        unless (exists $element_sequences{$data[0]}) { die "ERROR: Cannot find element $data[0]\n" } # a check to make sure the elemnent has been recorded
                        # check the orientation of the element sequence compared to the protein, report all sequences in the same
                        # orientation as the protein. Specifically test if the subject start position is lower than the end position or not
                        # to determine orientation.
                        if ($data[2] < $data[3]) {
                            $complete_elements_sequences{$data[0]} = $element_sequences{$data[0]};
                        }
                        else { # this element in on the other strand so change the orientation and reverse order of the positions on the location name
                            if ($data[0] =~ /(\S+):(\d+)-(\d+)/) {
                                my $updated_genomic_location = "$1:$3-$2";
                                $complete_elements_sequences{$updated_genomic_location} = rc($element_sequences{$data[0]});
                            }
                            else {
                                die "ERROR: Genomic location $data[0] could not be parsed $!";
                            }
                        }  
                    }
                }   

                # align the sequences 
                if (%complete_elements_sequences) { # only continue if complete elements were found
                    my $alignment_input_filename = File::Temp->new(UNLINK => 1);    # alignement program input and output files
                    my $alignment_output_filename = File::Temp->new(UNLINK => 1); 

                    # write sequences to the alignent input file, minus the TSDs
                    open (ALIINPUT, ">", $alignment_input_filename) or die "$!";
                    foreach my $header (keys %complete_elements_sequences) {

                        # modify the header to indicate if the TSDs are the same 
                        my $left_tir = substr($complete_elements_sequences{$header}, 0, $reviewed_tsdtirs{$cluster_name}[0]); 
                        my $right_tir = substr($complete_elements_sequences{$header}, length($complete_elements_sequences{$header}) - $reviewed_tsdtirs{$cluster_name}[0], $reviewed_tsdtirs{$cluster_name}[0]);                     
                      
                        my $element_sequence = substr($complete_elements_sequences{$header}, $reviewed_tsdtirs{$cluster_name}[0], length($complete_elements_sequences{$header})-(2*$reviewed_tsdtirs{$cluster_name}[0]));                       
                        print ALIINPUT ">$header\n$element_sequence\n";
                    }

                    close (ALIINPUT);

                    # run the alignment and convert result into a hash
                    `mafft --quiet $alignment_input_filename > $alignment_output_filename`;
                    if ($?) { die "ERROR executing mafft: error code $?\n"}  

                    # report results
                    my $alignment_file_output_name = "$cluster_name" . "_complete-elements-alignment.fa";
                    `cp $alignment_output_filename $ELEMENT_FOLDER/$cluster_name/$alignment_file_output_name`;
                    if ($?) { die "ERROR using cp: error code $?\n"}  

                    my $datestring = localtime(); 
                    print README "$datestring, File $alignment_file_output_name contains the alignment of nearly-complete elements\n";
                }
                else {
                    my $datestring = localtime(); 
                    print README "$datestring, in STEP 4, while TIRs were found no sequence matched the intial protein sequence\n";
                    print REJECT "$datestring\t$cluster_name\tSTEP 4\tNo matches to initial protein sequence with tblastn\n";
                    `mv $ELEMENT_FOLDER/$cluster_name $reject_folder_path`;
                    if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$cluster_name to $ANALYSIS_FOLDER: error code $?\n"}
                }
            }
            else {
                # my $datestring = localtime(); 
                # print README "$datestring, No genomic sequences with appropriately positioned TIRs were identified in STEP 4\n";
                # print REJECT "$datestring\t$cluster_name\tSTEP 4\tNo genomic sequences with appropriately positioned TIRs found\n";
                # `mv $ELEMENT_FOLDER/$cluster_name $reject_folder_path`;
                # if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$cluster_name to $ANALYSIS_FOLDER: error code $?\n"}
            }

            close README;
            print STDERR "done with $cluster_name\n";
        }
    }
}
close ANALYSIS;
close REJECT;


###### SUBSCRIPTS ##################
### These are scripts required for TE-discovery 

### Takes a fasta file name as input and return a hash with content
sub fastatohash {
    (my $filename) = @_;
    my %seqhash; #final hash with the seqhash
    my $current_header; # fasta header of the currently read sequence

    open (INPUT, $filename) or die "ERROR: cannot open input file $filename in fastatohash subroutine\n";
    my $line = <INPUT>;

    ## record the first header, ignoring everything after the first space in the header (to simplify recording later)
    if ($line =~ /^>(\S+)/) {
        $current_header = $1;
    }
    else {
        die "ERROR: file $filename does not appear to be a FASTA formatted file (or there's a space in the header after the > sign), this in the fastatohash subroutine\n";
    }

    while ($line = <INPUT>) {
        if ($line =~ /^>(\S+)/) {
            $current_header = $1;
            chomp $current_header; 
        }
        else {
            $line =~ s/\s//g; #remove all white spaces from the data line
            $seqhash{$current_header} .= $line;
        }
    } 

    return (%seqhash)
}

# Take a nucleotide sequence, location info and TSD length and return 1 if a TSD is found
sub gettsd {
	my ($seq, $loc1, $loc2, $type) = @_;
	### load the whole side sequences (to deal with gaps)
	my $seq_left_side = substr($seq, 0,  $loc1-1);
	my $seq_right_side = substr($seq, $loc2, -1);
	$seq_left_side =~ s/-//g; # remove gaps
	$seq_right_side =~ s/-//g;
    my $TSD_length;
    if ($type eq "TA") {
        $TSD_length = 2;
    }
    else {
        $TSD_length = $type;
    }

    my $c1 = cleanup(substr($seq_left_side, -$TSD_length, $TSD_length));
	my $c2 = cleanup(substr($seq_right_side, 0, $TSD_length));

    if ($type eq "TA") {
        if (($c1 eq "ta") and ($c2 eq "ta")) {
            return (1); # found a TA TSD
        }
        else {
            return (0);
        }
    }
    elsif (($type eq "2") and ($c1 and $c2)) { # check that both are equal and not 0, 0 means the sequence contained an "n" charcater
        if ($c1 eq "ta") { # this is to avoid duplication with the TA TSDs
            return (0); # this is a TA tsd, not approriate for this category
        }
        else {
            return (1); # found a 2 bp TSD that is not TA
        }
    }
    elsif (($c1 and $c2) and ($c1 eq $c2)) { # check that both are equal and not 0, 0 means the sequence contained an "n" charcater
            return (1); # found a TSD
    }
    else {
        return (0); # no TSD
    }
}

# given a sequence, locations of tirs, and maximum allowed mismatches between TIRs, returns the sequence of the tir if they are found, returns blanks otherwise
sub gettir {
	my ($seq, $loc1, $loc2, $min_tir_size, $max_number_mismatches) = @_;

    my $min_tir_size = 10; # anything smaller than this will not be reported as a TIR

    my $endfound = 0; #boolean 0 until the end of the TIR is found
	my $pos = 0; #current position in the sequence
	my $lastgoodbase = 0; #position of the last match of bases
	my $miss = 0; #number of non-matching sequences

	### load the sequence into memory and remove gaps
	my $sequence = substr($seq,$loc1-1,$loc2-$loc1+1); # DNA sequence of the whole element
    $sequence =~ s/-//g;
    ### get the ends into string Variables
    my $number_of_bp_to_scan = int((length $sequence)/2);
	my $s1 = substr ($sequence, 0, $number_of_bp_to_scan);
	my $s2 = substr ($sequence, -$number_of_bp_to_scan, $number_of_bp_to_scan);

    while (($pos <= (length $s1)) && ($endfound == 0)) {
        my $leftbase = substr($s1, $pos, 1); #base on the left end
        my $rightbase = substr($s2, -$pos -1, 1);

        # figure out if the bases match and are regular bases (not N)
        if (($leftbase eq (rc($rightbase))) and (("acgt" =~ /$leftbase/i) and ("acgt" =~ /$rightbase/i))) {
                $lastgoodbase = $pos;
        }
        else {
            $miss++;
        }

        #take stock if we need to stop
        if ($miss > $max_number_mismatches) {
            $endfound = 1;
        }
        else {
            $pos++;
        }
	}

    ### get the TIR sequences if a long enough tir has been found
    if ($min_tir_size <= $lastgoodbase) {
        my $tir1_sequence = substr($s1, 0, ($lastgoodbase+1));
        my $tir2_sequence = substr($s2, -($lastgoodbase+1), ($lastgoodbase+1));
        return ($tir1_sequence, $tir2_sequence);
    }
    else {
        return ("","");
    }
}

# Takes a string, converts it to lower case and checks if it contain an "n"
sub cleanup {
	my ($s) = @_;
    $s = lc($s);
    if ($s =~ /n/) {
        return (0); # this sequence contain at least one n character
    }
    else {
        return ($s);
    }
}

# Takes a folder name and standardizes the formating
sub fixdirname {
	my ($string) = @_;
	chomp $string;

	# if the name starts with ./ replace the period with the current directory
	if(substr($string, 0, 2) eq "./") {
		my $currdir = `pwd`;
		chomp $currdir;
		$currdir =~ s/ /\\ /g;
		$string = $currdir . "/" . substr($string, 2);
	}

	# replace spaces with '\ ' if they don't already have one
	$string =~ s/\\ /backslachandspace/g; # replace all existing backslash and spaces with long, unique, word
	$string =~ s/ /\\ /g; #replace all remaing spaces with backslash and space symbols
	$string =~ s/backslachandspace/\\ /g; #put the backslash and space symbols back
	$string =~ s/\(/\\(/g; #replace all open parentheses
	$string =~ s/\)/\\)/g; #replace all close parentheses

	#make sure the name does not end with a /
	if ((substr $string, -1) eq "/") {
		$string = substr($string, 0, -1);
	}
	return ($string);
}

#reverse complement
sub rc {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
}

# identify element sequences
sub identify_element_sequence {
    my ($chr_seq, $tir1, $tir2, $maximum_size, $chromosome_name, $tsd_length, $orientation) = @_; # takes as input the chromosome sequence, the TIR sequences,
                                                                                     # the maximum element size of the element, the chrosome name and orienation of the input sequence
    my @tir1; # location of all TIR1 sequence
    my @tir2; # location of all TIR2 sequence 
    my %nucleotide_sequences; # location of elements as key and nucleotide sequence as value, this is what is returned 
    # test that the orientation provided is known
    unless (($orientation eq "+") or  ($orientation eq "-")) { die "ERROR: Orientation $orientation is not known in subroutine identify_element_sequences\n$tir1, $tir2, $maximum_size, $chromosome_name, $orientation\n"}

    # if the search says to search for the reverse strand, change the tir sequences
    if ($orientation eq "-") {
        $tir1 = rc($tir2);
        $tir2 = rc($tir1);
    }
    while ($chr_seq =~ m/$tir1/g) {
        push @tir1, pos($chr_seq) - length($tir1); # adjust the left position to compensate for where the loc is found
    }
    while ($chr_seq =~ m/$tir2/g) {
        push @tir2, pos($chr_seq);
    }
    # identify any pairs that are in the correct orientation and not too far from one another
    for (my $i=0; $i<scalar @tir1; $i++) {
        for (my $j=0; $j<scalar @tir2; $j++) {
            my $size = $tir2[$j] - $tir1[$i] + $tsd_length + $tsd_length;
            if (($tir1[$i] < $tir2[$j]) and ($size < $maximum_size)){
                my $b1; # boundary of the element on the chromosome
                my $b2; # boundary of the element on the chromosome
                if ($orientation eq "+") {
                    $b1 = $tir1[$i]+1-$tsd_length; 
                    $b2 = $tir2[$j]+$tsd_length;
                }
                else {
                    $b2 = $tir1[$i]+1-$tsd_length;
                    $b1 = $tir2[$j]+$tsd_length;
                }

                # write the output and check if the TSDs are the same and update the location name if they are
                my $location_name = "$chromosome_name:$b1-$b2"; 
                my $full_sequence = substr($chr_seq, $tir1[$i]-$tsd_length, $size);
                my $tsd1 = substr($full_sequence, 0, $tsd_length);
                my $tsd2 = substr($full_sequence, length($full_sequence)-$tsd_length, $tsd_length);
                if ($tsd1 eq $tsd2) {
                    $location_name .= "-identicalTSDs($tsd1)";
                }
                $nucleotide_sequences{$location_name} = $full_sequence;
            }
        }
    }
    return (%nucleotide_sequences);
}

# in-depth comparion between possible tirs and what's already been seen
# returns 2 if it's an exact match, 1 for a partial match, and 0 for no match. Postitive numbers 
# if the matches are in the orientation, negative for reverse orientation
sub compare_tirs {
    my($t1, $t2, %seen) = @_; # tirs to explore and what's been seen
    foreach my $elementname (keys %seen) {
        my $s1_match=0; 
        my $s2_match=0; 

        # test in the same orientation                
        $s1_match = testmatch($t1,$seen{$elementname}[1],1); # the inputs are the two sequences to compare, and 1 to indicate it's the left tir   
        $s2_match = testmatch($t2,$seen{$elementname}[2],2); # the inputs are the two sequences to compare, and 2 to indicate it's the right tir 
        if ($s1_match and $s2_match) { # match found
            if (($t1 eq $seen{$elementname}[1]) and ($t2 eq $seen{$elementname}[2])) { # true if the TIRs are exactly the same
                return(2, $seen{$elementname}[1], $seen{$elementname}[2], $seen{$elementname}[0]);

            }
            else { # the match is not exact
                return(1, $seen{$elementname}[1], $seen{$elementname}[2], $seen{$elementname}[0]);
            }
        }

       # test in the reverse orientation
       $s1_match = testmatch(rc($t2),$seen{$elementname}[1],1); 
       $s2_match = testmatch(rc($t1),$seen{$elementname}[2],2); 
       if ($s1_match and $s2_match) {
            if ((rc($t1) eq $seen{$elementname}[2]) and (rc($t2) eq $seen{$elementname}[1])) { # true if the TIRs are exactly the same
                return(-2, $seen{$elementname}[1], $seen{$elementname}[2], $seen{$elementname}[0]);
            }
            else {
                return(-1, $seen{$elementname}[1], $seen{$elementname}[2], $seen{$elementname}[0]);
            }
       }
    }
    return (0,0,0);

    # compare a pair of TIRs
    sub testmatch {
        my ($s1, $s2,$side) = @_; # the two input sequences and an indication of which side comparing
        if ((length $s1) < (length $s2)) {
            if ($s2 =~ /$s1/) {
                if (($-[0] == 0) and ($side == 1)) { # test if the match starts at the first position of the left TIR
                    return(1);
                }
                elsif (($+[0] == length($s2)) and ($side == 2)) { # test if the match starts at the first position of the right TIR
                    return(1)
                }
            }
        }
        else {
            if ($s1 =~ /$s2/) {
                if (($-[0] == 0) and ($side == 1)) { # test if the match starts at the first position of the left TIR { # test if the match starts at the first position of the TIR
                    return(1)
                }
                elsif (($+[0] == length($s1)) and ($side == 2)) { # test if the match starts at the first position of the right TIR
                    return(1)
                }
            }
        }
        return(0);
    }  
}

sub average_tir_number_and_length {
    my ($loc1, $loc2, $min_tir_size, $max_number_mismatches, %seq) = @_; # get the bounds, size, and sequences
    my $number_of_sequences; # total sequences in this alignment
    my $total_TIR_length; # sum of all the TIR length, used to calculate the average
    my $TIR_number; # number of sequences with TIRs
    my $average_TIR_length;

    foreach my $sequence_name (keys %seq) {
        unless ($sequence_name =~ /consensus/) {
            my ($TIR1, $TIR2) = gettir($seq{$sequence_name}, $loc1, $loc2, $min_tir_size, $max_number_mismatches);
            if ($TIR1) {
                $total_TIR_length += length $TIR1;
                $TIR_number += 1;
            }
        }
    }
    
    if ($TIR_number) { # if TIRs have been found
        $average_TIR_length = int($total_TIR_length / $TIR_number);
    }
    else {
        $TIR_number = "No TIRs found";
        $average_TIR_length = "N/A";
    }

    return ($TIR_number, $average_TIR_length);
}