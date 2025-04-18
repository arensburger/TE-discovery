# TE Discovery Pipeline - Main Script
# Author: Peter Arensburger
# Date: March 2025
# Description: This script runs the full TE discovery pipeline using input protein sequences and genome data.

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use FindBin::libs;  # this sets up that the directory lib will have all modules necessary to run the program 
use TEdiscovery;    # these are the subcripts necessary to run the pipeline
use List::UtilsBy qw(max_by);
use List::Util qw(max);
use Scalar::Util qw(looks_like_number);

### CONSTANTS for all steps
my $NUM_THREADS = `nproc --all`;
if ($?) { warn "WARNING could not determine the number of cores automatically, defaulting to 8\n"; $NUM_THREADS=8}
chomp $NUM_THREADS;

### INPUTs from command line, top level variables are in uppercase
my $INPUT_PROTEIN_SEQUENCES; # fasta formated file with input protein sequences
my $TBLASTN_FILE; # name of the file containing the out put of the tblastn
my $INPUT_GENOME; # fasta formated file with genome that input proteins
my $ANALYSIS_NAME; # name to give to this analysis, can be any string
my $ANALYSIS_FILES_OUTPUT_DIR; # directory where analsysis output files of the analysis are stored
my $ELEMENT_FOLDER; # directory where individual folders for each element are stored
my $START_STEP = 0; # analysis step to start at, by default is set to zero
my $END_STEP = 100000; # analysis step to end at, by deault set to 1000000 (hopefully fewer than those number of steps)

### CHECK INPUTS Read and check that the inputs have been provided
GetOptions(
	'p:s'   => \$INPUT_PROTEIN_SEQUENCES,
	't:s'   => \$TBLASTN_FILE,
	'g:s'   => \$INPUT_GENOME,
    'n:s'   => \$ANALYSIS_NAME,
    'a:i'   => \$START_STEP,
    'b:i'   => \$END_STEP,
);

## CHECK INPUTS Validate required input files
## There are two possible combinations of inputs
unless (($INPUT_PROTEIN_SEQUENCES and $INPUT_GENOME and $ANALYSIS_NAME) or ($TBLASTN_FILE and $INPUT_GENOME and $ANALYSIS_NAME)) {
	die "Usage, two options:
    Option 1: perl TEdiscovery.pl <-p protein fasta file REQUIRED> <-g fasta file of genome processed with makeblastdb as shown below REQUIRED> <-n name for this analysis REQUIRED> <-a starting analysis step OPTIONAL> <-b end analysis step OPTIONAL>\n
    Option 2: perl TEdiscovery.pl <-t output of tblastn REQUIRED> <-g fasta file of genome processed with makeblastdb as shown below REQUIRED> <-n name for this analysis REQUIRED> <-a starting analysis step OPTIONAL> <-b end analysis step OPTIONAL>
    
    Commands to generate required files for both options:
    1) makeblastdb -in <genome fasta file name> -dbtype nucl
    2) samtools faidx <genome fasta file name>
       awk \'{OFS=\"\\t\"; print \$1,\$2}\' < <genome fasta file name>.fai > <genome fasta file name>.length

    Command to generate required file for option 2
    tblastn -query <protein fasta file> -db <genome fasta file name> -outfmt \"6 qseqid sseqid sstart send pident length qlen\" -out tblastn.o -num_threads <number of threads>
    (the output file of this command should be provided as the -t parameter)\n";
}

## VARIABLES used by more than one step in the pipeline
my $ANALYSIS_FILES_OUTPUT_DIR="./$ANALYSIS_NAME-analysis-files"; # directory to store output files of current analysis (no slash at the end), these can be destroyed when analysis is finished

## CHECK INPUTS Create output directory for analysis files if necessary
print "Preliminary steps...\n";
if (-d $ANALYSIS_FILES_OUTPUT_DIR) {
    die "Directory $ANALYSIS_FILES_OUTPUT_DIR already exists, the script needs to create a new empty directory\n";
}
else {
    print "\tCreating directory $ANALYSIS_FILES_OUTPUT_DIR for storing files generated during the analysis\n";
    mkdir( $ANALYSIS_FILES_OUTPUT_DIR ) or die "Couldn't create $ANALYSIS_FILES_OUTPUT_DIR directory, $!";
}

## CHECK INPUTS Create output directory for individual elements if necessary
my $ELEMENT_FOLDER="./$ANALYSIS_NAME-elements"; # directory where the analysis of individual element will be stored
if (-d $ELEMENT_FOLDER) {
    die "Directory $ELEMENT_FOLDER already exists, the script needs to create a new empty directory\n";
}
else {
    print "\tCreating directory $ELEMENT_FOLDER that will have subdirectories for individual elements\n";
    mkdir( $ELEMENT_FOLDER ) or die "Couldn't create $ELEMENT_FOLDER directory, $!";
}

## Create files to store analysis parameters, and to store rejected sequences. 
## Search for a unique names, with the same index for both files
my $analysis_parameters_file_name = "$ANALYSIS_FILES_OUTPUT_DIR/Analysis_parameters.txt"; # file to record parameters
my $rejection_file_name = "$ANALYSIS_FILES_OUTPUT_DIR/Rejected_sequences.txt"; # file to store rejected sequences, and why
if ((-f $analysis_parameters_file_name) or (-f $rejection_file_name)){
   die "ERROR: Either or both files $analysis_parameters_file_name and/or $rejection_file_name already exists, this should not happen"
}
else {
    open (ANALYSIS,'>', $analysis_parameters_file_name) or die "ERROR: cannot open file $analysis_parameters_file_name\n";
    open (REJECT, '>', $rejection_file_name) or die "ERROR, cannot create output file $rejection_file_name\n";
}
print ANALYSIS "Analysis name: $ANALYSIS_NAME\n";
my $datestring = localtime();
print ANALYSIS "Date and time: $datestring\n";
print ANALYSIS "Input file: $INPUT_PROTEIN_SEQUENCES\n";
print ANALYSIS "Genome: $INPUT_GENOME\n\n";
print ANALYSIS "Parameters, set in the script:\n";

my $BLAST_OUTPUT_FILE_NAME = "$ANALYSIS_FILES_OUTPUT_DIR/tblastn.o"; # default name and location unless a file is provided, default file name has index added
                                                                                # also, could not define this earlier because it depends on the ultimate analysis name

### PIPELINE STEP 1 identify proteins that match the genome with parameters specified above under
###     The output is a list of proteins for further analysis recorded in the file $output_file_name
### CONSTANTS applicable to this step only (also record these in the file)
my $GENOME_IDENTITY = 80; # IDENTIFYING PROTEINS, per protein, minimum percent identity between protein and genome
print ANALYSIS "STEP1: GENOME_IDENTITY = $GENOME_IDENTITY\n";
my $COVERAGE_RATIO = 0.5; # IDENTIFYING PROTEINS, per protein, minimum ratio of (blast match length) / (query length)
print ANALYSIS "STEP1: COVERAGE_RATIO = $COVERAGE_RATIO\n";
my $COPY_NUMBER = 2; # IDENTIFYING PROTEINS, minimum number of copies that hit different parts of the genome 
print ANALYSIS "STEP1: COPY_NUMBER = $COPY_NUMBER\n";
my $MIN_DISTANCE = 10000;   # IDENTIFYING PROTEINS, if two elements are on the same chromosome, how far they have to be, to be considered different elements
                            # NOTE: the minimum distance should bigger than the $BLAST_EXTEND variable, to avoid having the same element recorded twice        
print ANALYSIS "STEP1: MIN_DISTANCE = $MIN_DISTANCE\n";

my $step_number = 1;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

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
        print ANALYSIS "STEP1: tblastn output was provided in file $TBLASTN_FILE\n";
    }
    else {
        `tblastn -query $INPUT_PROTEIN_SEQUENCES -db $INPUT_GENOME -outfmt "6 qseqid sseqid sstart send pident length qlen" -out $BLAST_OUTPUT_FILE_NAME -num_threads $NUM_THREADS`;
        if ($?) { die "ERROR executing tblastn, stopping analysis (hint: was the genome formated with makeblastdb?): error code $?\n"}
        print ANALYSIS "STEP1: tblastn was run by the scritpt the output is in file $BLAST_OUTPUT_FILE_NAME\n";
    }

    ## Inspired by the Goubert et al. protocol, filter elements that 1) have >= 80% identity to genome, 2) have 50% length of the query, 3) are found at multiple locations
    my %candidate_protein; # hash with protein name as key and string with chromosome and middle location of element on that chromsome
    open (INPUT, "$BLAST_OUTPUT_FILE_NAME") or die "ERROR: Cannot open file $BLAST_OUTPUT_FILE_NAME\n";
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
        }   
        else {
            $rejected_ids{$data[0]} = "STEP $step_number\tError code (genome identity/coverage/minimum distance): $gi $cr $md"
        }   
    }
    close INPUT;

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
            print REJECT "$datestring\t$prot_name\tSTEP $step_number\tError code (BLAST minimum copy number/observed copy number) $COPY_NUMBER $protein_ids{$prot_name}\n";
            delete $protein_ids{$prot_name};
        }
    }
    close OUTPUT;

    # Create individual directories for each element
    my $i=0; # counts the number of output lines, to check if it's zero
    foreach my $prot_name (keys %protein_ids) {
        if (-d "$ELEMENT_FOLDER/$prot_name") {
            print "\tWARNING: Element folder $ELEMENT_FOLDER/$prot_name already exists (this should not normally happen), writing files to this folder\n";
        }
        else {
            mkdir( "$ELEMENT_FOLDER/$prot_name" ) or die "Couldn't create $ELEMENT_FOLDER/$prot_name directory, $!";
        }
        $i++;
    }

    if ($i) {
        print "Finished STEP $step_number, identified $i candidates for further analysis\n";
    }
    else {
        warn "WARNING: STEP $step_number did not result in any identified candiates, no output produced\n";
    }
}

### PIPELINE STEP 2 
### For each element that had enough approved blast hits, look for TSD-TIR <---> TIR-TSD combinations 
### CONSTANTS applicable only for STEP 2
my $BLAST_EXTEND = 2000; # IDENTIFYING PROTEINS, number of bp to extend on each side of the blast hit
print ANALYSIS "STEP2: BLAST_EXTEND = $BLAST_EXTEND\n";
my $MAX_SEQUENCE_NUMBER = 100; # ALIGNING SEQUENCES maximum number of sequences to consider, to save time
print ANALYSIS "STEP2: MAX_SEQUENCE_NUMBER = $MAX_SEQUENCE_NUMBER\n";
my $GAP_THRESHOLD=0.75; # REMOVING GAPS FROM ALIGNMENT, if an alignment position has this proportion or more of gaps, then remove it from the multiple sequence alignment
print ANALYSIS "STEP2: GAP_THRESHOLD = $GAP_THRESHOLD\n";
my $CONSLEVEL=0.60; # MAKING CONSENSUS OF SEQUENCES sequence consensus level for consensus
print ANALYSIS "STEP2: CONSLEVEL = $CONSLEVEL\n";
my $WINDOW_SIZE = 15; # MAKING CONSENSUS OF SEQUENCES size of the window looking for stretches of N's
print ANALYSIS "STEP2: WINDOW_SIZE = $WINDOW_SIZE\n";
my $SCAN_PROP = 0.5; # MAKING CONSENSUS OF SEQUENCES minimum proportion of side scan that has to be N's to be a real N full edge
print ANALYSIS "STEP2: SCAN_PROP = $SCAN_PROP\n";
my $MAX_WIN_N = 2; # MAKING CONSENSUS OF SEQUENCES maximum number of N's in the first window where the transition from N to non-N is
print ANALYSIS "STEP2: MAX_WIN_N = $MAX_WIN_N\n";
my $EDGE_TEST_PROPORTION = 0.05; # TESTING CONSENSUS SEQUENCES how far from the edge of the consensus do the non-gap positions have to start
print ANALYSIS "STEP2: EDGE_TEST_PROPORTION = $EDGE_TEST_PROPORTION\n";
my $MIN_TIR_SIZE = 10; # IDENTIFYING TIR-TSDS smallest allowable size for the TIR
print ANALYSIS "STEP2: MIN_TIR_SIZE = $MIN_TIR_SIZE\n";
my $TIR_MISMATCHES = 2; # IDENTIFYING TIR-TSDS maximum number of mismatches allowed between two TIRs
print ANALYSIS "STEP2: TIR_MISMATCHES = $TIR_MISMATCHES\n";
my $TIR_PROP_CUTOFF = 0.15; # IDENTIFYING TIR-TSDS proportion of elements with TIRs at which positions are reported
print ANALYSIS "STEP2: TIR_PROP_CUTOFF = $TIR_PROP_CUTOFF\n";
my $MIN_PROPORTION_SEQ_WITH_TIR=0.25; #IDENTIFYING TIRs minimum proportion of total elements for a sequence that must contain proper TIRs to be considered a candidate
print ANALYSIS "STEP2: MIN_PROPORTION_SEQ_WITH_TIR = $MIN_PROPORTION_SEQ_WITH_TIR\n";
my $MAX_TIR_PROPORTION=0.75; #IDENTIFYING TIRs how close to the maximum number of tirs do you have to be to qualify as a top TIR
print ANALYSIS "STEP2: MAX_TIR_PROPORTION = $MAX_TIR_PROPORTION\n";
my $MAX_END_PROPORTION=0.75; #IDENTIFYING TIRs how close to maximum proportion of sequences with identical start and stop of tir sequences you can be to a top number
print ANALYSIS "STEP2: MAX_END_PROPORTION = $MAX_END_PROPORTION\n";
my $MAX_TSD_PROPORTION=0.5; #IDENTIFYING TIRs how close to maximum number of TSDs to qualify as a top TSD sequence
my %EXAMINE_CODES=("1111" => 1, "1101" => 2); # success codes to examine as key and priority as value
print ANALYSIS "STEP2: %EXAMINE_CODES=(\"1111\" => 1, \"1101\" => 2) # TIR-TSD success codes that will allow an element to be examine further as key and priority as value\n";
print ANALYSIS "STEP2: MAX_TSD_PROPORTION = $MAX_TSD_PROPORTION\n";
print ANALYSIS "STEP2: TIR-TSD success code, first number: Is the proportion of sequences with TIRs higher than MIN_PROPORTION_SEQ_WITH_TIR?\n";
print ANALYSIS "STEP2: TIR-TSD sode, second number: Does this candidate have one of the highest number of TIRs?\n";
print ANALYSIS "STEP2: TIR-TSD sode, third number: Do the TIRs tend to start and end with the same bases?\n";
print ANALYSIS "STEP2: TIR-TSD sode, first number: Does this candidate have one of the highest number of TSDs?\n";

my $step_number = 2;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## VARIABLES, variable for this step
    my @elements; # name of all the elements that will be anlaysed in this step

    ## STEP 2.1
    ## figure out the elements that we're working with 
    ## (for the sake of being modular, redoing this instead of just taking the data from the previous step,
    ## this also ensures that a folder has been created for each element).

    opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
    while (readdir $dh) {
        unless ($_ =~ /^\./) { # prevents reading invisible files or the . and .. files
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
        print "\tProcessing $element_name, number $i of $size\n";
        $i++;

        # create or open the README file for this element
        open (README, ">$ELEMENT_FOLDER/$element_name/README.txt") or die "ERROR: Could not open or create README file $ELEMENT_FOLDER/$element_name/README.txt\n";

        # STEP 2.2.1
        # for each sequence of this element extend it by $BLAST_EXTEND bps on both sides of the sequence
        my @blastlines = (); # holds all the relevant blast lines for this element
        my $j; # counter of the number of blast lines for this element
        open (INPUT, $BLAST_OUTPUT_FILE_NAME) or die "ERROR: cannot open file blast file output $BLAST_OUTPUT_FILE_NAME\n";
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

        # checking that the genome length file is present
        unless (-f "$INPUT_GENOME.length") {
            die "ERROR: Cannot file find file $INPUT_GENOME.length needed for bedtools, generate one using\n(assuming it's been formated for BLAST with makeblastdb)\n\nsamtools faidx \$genome\nawk \'{OFS=\"\\t\"; print \$1,\$2}\' < \$genome.fai > \$genome.length\n";
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

         my $extended_fasta_name = File::Temp->new(UNLINK => 1, SUFFIX => '.fa' ); # name fo the file with extended fasta sequences
        `bedtools getfasta -fi $INPUT_GENOME -fo $extended_fasta_name -bed $slopfile -s`;
        if ($?) { die "ERROR executing bedtools: error code $?\n"}

        # STEP 2.2.2
        # align the sequences
        my $aligned_sequences = File::Temp->new(UNLINK => 1, SUFFIX => '.maf' ); 
     
       `mafft --quiet --thread -1 $extended_fasta_name > $aligned_sequences`;
        if ($?) { die "Error executing mafft, error code $?\n"}

        # STEP 2.2.3
        # remove positions from the alignment that have more gaps than the threshold given $GAP_THRESHOLD   
        # a record of the deletions is made so the original alignment can be recreated later

        my %cp; # cp = current position, sequence name as key and current position as value at [0] and orientation at [1]
                # [2] is boolean, value 0 until a non-gap has been reached
        my %seq; # sequence name as key and string with mafft alignment as value
        my %seqrmg; # sequences with gaps removed
        my %pos; # holds the sequence name as key array with position of every nucleotide as value
        my $alilen; # length in bp of the multiple sequence alignment
        my $ali_trimmed_length; # length in bp of the multiple squence alignment after it's been trimmed

        # populate the %cp hash with information about the multiple sequence alignment
        my %seq = fastatohash($aligned_sequences); 
        foreach my $name (keys %seq) { 
            if ($name =~ /^(\S+):(\d+)-(\d+)\((.)\)/) {
                my $chr = $1;
                my $b1 = $2;
                my $b2 = $3;
                my $ori = $4;
                if ($ori eq "+") {
                    $cp{$name}[0] = $b1;
                    $cp{$name}[1] = 1;
                    $cp{$name}[2] = 0;
                }
                elsif ($ori eq "-") {
                    $cp{$name}[0] = $b2;
                    $cp{$name}[1] = -1;
                    $cp{$name}[2] = 0;
                }
                else {
                    die "ERROR: Reading sequence alignment for element $element_name, orientation not recognized on this element $name\n";
                }
            }
            else {
                die "ERROR: Reading sequence alignment for element $element_name, fasta title not formated as expected\n$name";
            }
            $alilen = length $seq{$name}; # this is the length of the sequence alignment prior to trimming
        }

        # go through each position of the alignment and decide if the position needs to be removed
        for (my $i=0; $i<$alilen; $i++) { # go through positions of the alignment
            my $numgap; # number of gaps in this column
            foreach my $name (keys %seq) { # go through each element of the mafft alignment
                if (substr($seq{$name},$i,1) eq "-") {
                    $numgap++;
                }
                else { # if this position does not have a gap then update %cp
                    if ($cp{$name}[3] == 0) { # if this is the first time the position has a non-gapped position
                        $cp{$name}[3] = 1;
                    }
                    else { # not the first encounter and position should be updated
                        $cp{$name}[0] += $cp{$name}[1]; # position will be increased or decreased based on orientation
                    }
                }
            }
            unless ($numgap/(keys %seq) >= $GAP_THRESHOLD) { # if this column is kept, then copy it to %seqrmg and update %pos
                foreach my $name (keys %seq) { # go through each element of the mafft alignment
                    $seqrmg{$name} .= substr($seq{$name},$i,1);
                    push @{ $pos{$name} },$cp{$name}[0];
                }
                $ali_trimmed_length++;
            }
        }

        # printing the .maf later, after the consensus has been created

        # output the file that will allow restoring the deleted positions, called .alipos
        open (OUTPUT, '>', "$ELEMENT_FOLDER/$element_name/$element_name.alipos") or die "Error: cannot create file $element_name.alipos, $!\n";
        for my $name (keys %pos) {
            for my $j ( 0 .. $#{ $pos{$name}}) {
                my $location = $j+1;
                print OUTPUT "$name\t$location\t$pos{$name}[$j]\n";
            }
	    }
        my $datestring = localtime(); # update the README file of new file created
        print README "$datestring, File $element_name.alipos contains the information about positions removed in the alignment file\n";
	    close OUTPUT;

        # STEP 2.2.4
        # Create a consensus sequence for this element
        my $conseq; # initial consensus sequence (not trimmed on the ends)
        my @prop_cons; # for each nucleotide position of the consensus sequence holds the proportion of sequences that have the most abundant nucleotide 
        my @prop_complete_total; # for each nucleotide position the proportion of the complete total that have the consensus

        # loop though all positions and make a consensus sequence
        for (my $i=0; $i < $ali_trimmed_length; $i++) {
            my %abundance; # holds the nucleotide identify as key and frequency as value
            my $total; # total number of relevant squences for this column
            foreach my $key (keys %seqrmg) { # loop through all sequences and record any identified nucleotides
                my $character = lc(substr($seqrmg{$key}, $i, 1));
                if (($character eq "a") or ($character eq "c") or ($character eq "g") or ($character eq "t")) {
                    $abundance{$character} += 1; # add the sequences
                }
                $total++; # count the number of sequences
            }

            # add the highest abundance nucleotide if above consensus level threshold
            my $highest = max_by { $abundance{$_} } keys %abundance; # from https://perlmaven.com/highest-hash-value           
            if (($abundance{$highest}/$total) >= $CONSLEVEL) {
                $conseq .= $highest;
            }
            else {
                $conseq .= "N";
            }
        }

        # STEP 2.2.5
        # trim the ends of the consensus that have low agreement on a single sequence
        my $ltrans=0; # position of the transition from "not the element" to the "element" on the left side of the alignment
        my $rtrans=length $conseq; # position of the transition from "not the element" to the "element" on the right side of the alignment
        my $trimmed_conseq; # consensus sequence trimmed of side sequences, remains blank if no transition was found

        # record the number of N's in each window
        my @Nnum; # number of N's (or non-known bases) in the consensus in windows of size $WINDOW_SIZE
        for (my $i=0; $i <= ((length $conseq)-$WINDOW_SIZE); $i++) {
            my $winseq = substr ($conseq, $i, $WINDOW_SIZE);
            my $nBases = $winseq =~ tr/ACGTacgt//; # $nBases holds number known bases in the current window
            push @Nnum, (length $winseq) - $nBases; # this add the number of non-bases to the array @Nnum for current window
        }

        # find where the transtion from high number of N's in the consenus to low number occurs
        # finding the transition on the left side
        my $i=0;
        while (($ltrans==0) and ($i < ((length $conseq)-$WINDOW_SIZE + 1))) {
            if ($Nnum[$i] <= $MAX_WIN_N) { # find the first window from the left that has $MAX_WIN_N propotion of N
                                            # this is where the transition is
                my $s = substr($conseq, $i, $WINDOW_SIZE);
                unless ($s =~ /^N/i) { # prevent the transtion occuring at an N
                    $ltrans = $i+1;
                }
            }
            $i++;
        }

        # finding the transition on the right side
        my $i=(length $conseq)-$WINDOW_SIZE;
        while (($rtrans==(length $conseq)) and ($i >= 0)) {
            if ($Nnum[$i] <= $MAX_WIN_N) {
                my $s = substr($conseq, $i, $WINDOW_SIZE);
                unless ($s =~ /^N/i) { # prevent the transtion occuring at an N
                    $rtrans = $i+$WINDOW_SIZE;
                }
            }
            $i--;
	    }

        # if both a left and right transition are found, then create the trimmed consensus sequence and add it to the alignment file
        if (($ltrans > 0) and ($rtrans < (length $conseq))) { 
            for (my $i=0; $i<$ltrans-1; $i++) {
                $trimmed_conseq .= "-";
            }
            for (my $i=($ltrans-1); $i<($rtrans); $i++) {
                $trimmed_conseq .= substr($conseq, $i, 1);
            }
            for (my $i=($rtrans); $i<(length $conseq); $i++) {
                $trimmed_conseq .= "-";
            } 
        }

        # output the alignment file with positions removed, called .maf add the consensus if there is one
        open (OUTPUT, '>', "$ELEMENT_FOLDER/$element_name/$element_name.maf") or die "ERROR: Cannot create file $element_name.maf, $!\n";
        if (length($trimmed_conseq)) {
            print OUTPUT ">consensus-$ltrans-$rtrans\n";
            print OUTPUT "$trimmed_conseq\n"
        }
        foreach my $name (keys %seqrmg) {
            print OUTPUT ">$name\n$seqrmg{$name}\n";
            $ali_trimmed_length = length $seqrmg{$name}; # length of the alignment after trimming
        }
        my $datestring = localtime(); # update the README file of new file created
        print README "$datestring, File $element_name.maf is the alignment of sequences with positions containing $GAP_THRESHOLD proportion of gaps removed\n";
        close OUTPUT;

        # Test to see if this consensus sequence passes the "edge test". The trimmed consensus should not be too close 
        # to the edge of the alignment (this would suggest that the whole region, not just an element inside it, is conserved)
        my $edge_test; # boolean, will be set to 1 if the consensus is not too close to the edge
        if (($ltrans > ($EDGE_TEST_PROPORTION * (length $trimmed_conseq))) and ($rtrans < ((1 - $EDGE_TEST_PROPORTION) * (length $trimmed_conseq)))) {
            $edge_test = 1;
        }

        # update the README file
        my $datestring = localtime(); 
        if ($edge_test and $trimmed_conseq) {
            print README "$datestring, a trimmed consensus sequence was created that passed the \"edge test\"\n";
        }
        elsif ($trimmed_conseq) { # get here if a consenus was found but it failed the edge test
            print README "$datestring, a trimmed consensus sequence was created but it failed the \"edge test\", stopping the analysis here\n";
            print REJECT "$datestring\t$element_name\tSTEP 2\tfailed the EDGE TEST\n";
            `mv $ELEMENT_FOLDER/$element_name $ANALYSIS_FILES_OUTPUT_DIR`;
            if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FILES_OUTPUT_DIR: error code $?\n"}
        }
        else { # get here if a no abrupt changes in the consensus sequences was detected
            print README "$datestring, the consensus sequence showed no transitions into an element, stopping the analysis here\n";
            print REJECT "$datestring\t$element_name\tSTEP 2\tNo transitions to element observed in alignment file\n";
            `mv $ELEMENT_FOLDER/$element_name $ANALYSIS_FILES_OUTPUT_DIR`;
            if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FILES_OUTPUT_DIR: error code $?\n"}
        }

        # STEP 2.2.6 
        # identify the TIR and TSD locations around the edges of the consensus sequence
        if ($trimmed_conseq and $edge_test) { # only continue if consensus was created and edge test was passed

            # Variables relevant to identifying TIRs and TSDs
            my $range = 5; # how many bp to search around for tirs
            my $max_TIR_number; # highest number of TIRs observed for one pair of start and end positions
            my $max_proportion_first_last_bases; # highest number of locations that start and end with the same bases
            my $max_TSD_number; # highest number of intact TSDs
            my @tsd_tir_combinations; # array of possible locations along with various information

            for (my $i=$ltrans-$range; $i<=$ltrans+$range; $i++) {
	            for (my $j=$rtrans-$range; $j<=$rtrans+$range; $j++) {
		            my $number_of_tirs_found=0; # nubmer of sequences that match the TIR criteria
                    my %tir_first_and_last_bases; # first and last set of bases of tir as key and abundance as value
                    my %tsds_found; # keys is TSD type "TA", "2", ... "10" and key is number of TSDs found
                    foreach my $sequence_name (keys %seqrmg) {
                        my ($tir1_sequence, $tir2_sequence) = gettir($seqrmg{$sequence_name}, $i, $j, $MIN_TIR_SIZE, $TIR_MISMATCHES); # figure if this sequences has a tir at these positions and if so, report first and last nucleotide
                        if ($tir1_sequence) {
                            $number_of_tirs_found += 1;
                            my $bases = substr($tir1_sequence, 0, 3) . substr($tir2_sequence, -3, 3); # recording the first and last 3 bases 
                            $tir_first_and_last_bases{$bases} += 1;
                        }
                        $tsds_found{"TA"} += gettsd($seqrmg{$sequence_name}, $i, $j, "TA");
                        $tsds_found{2} += gettsd($seqrmg{$sequence_name}, $i, $j, 2);
                        $tsds_found{3} += gettsd($seqrmg{$sequence_name}, $i, $j, 3);
                        $tsds_found{4} += gettsd($seqrmg{$sequence_name}, $i, $j, 4);
                        $tsds_found{5} += gettsd($seqrmg{$sequence_name}, $i, $j, 5);
                        $tsds_found{6} += gettsd($seqrmg{$sequence_name}, $i, $j, 6);
                        $tsds_found{7} += gettsd($seqrmg{$sequence_name}, $i, $j, 7);
                        $tsds_found{8} += gettsd($seqrmg{$sequence_name}, $i, $j, 8);
                        $tsds_found{9} += gettsd($seqrmg{$sequence_name}, $i, $j, 9);
                        $tsds_found{10} += gettsd($seqrmg{$sequence_name}, $i, $j,10);
                    }

                    # figure out the most abundant begining and end of the TIR 
                    my $most_abundant_tir_seq=0;
                    foreach my $name (sort { $tir_first_and_last_bases{$a} <=> $tir_first_and_last_bases{$b} } keys %tir_first_and_last_bases) {
                        $most_abundant_tir_seq = $tir_first_and_last_bases{$name}/$number_of_tirs_found;
                    }

                    # record all the data for this compbination of $i and $j and update maximum values
                    $max_TIR_number = max($max_TIR_number, $number_of_tirs_found);
                    $max_proportion_first_last_bases = max ($max_proportion_first_last_bases, $most_abundant_tir_seq);
                    $max_TSD_number = max ($max_TSD_number, $tsds_found{"TA"}, $tsds_found{2}, $tsds_found{3}, $tsds_found{4}, $tsds_found{5}, $tsds_found{6}, $tsds_found{7}, $tsds_found{8}, $tsds_found{9}, $tsds_found{10});
                    push @tsd_tir_combinations, "$i\t$j\t$number_of_tirs_found\t$most_abundant_tir_seq\t$tsds_found{\"TA\"}\t$tsds_found{2}\t$tsds_found{3}\t$tsds_found{4}\t$tsds_found{5}\t$tsds_found{6}\t$tsds_found{7}\t$tsds_found{8}\t$tsds_found{9}\t$tsds_found{10}";
                }
            }

            # go through the element and identify those candidate locations that pass the tests for TIR-TSD combinatations
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
                if (($number_of_tirs/(keys %seqrmg)) >= $MIN_PROPORTION_SEQ_WITH_TIR) {
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
                `mv $ELEMENT_FOLDER/$element_name $ANALYSIS_FILES_OUTPUT_DIR`;
                if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FILES_OUTPUT_DIR: error code $?\n"}
            }
        }     
        close README;
    }    
}

### PIPELINE STEP 3 
### Present the elements to the user for manual review
### CONSTANTS applicable only for STEP 3
my $TIR_bp = 30; # how many bp to display on the TIR side
my $temp_alignment_file = "/tmp/ali.fa"; # tried using perl temporary file system, but aliview will not open those

my $step_number = 3;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    my $pkey; # pressed key, used for input from user

    ## make a list of elements to analyze, put those element into %files along with relevant file names
    my %files; # holds the element name as key and path to the .tirtsd and the .maf files as array of values
    opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
    while (readdir $dh) {
        unless ($_ =~ /^\./) { # prevents reading invisible files or the . and .. files
            my $specific_element_folder = $ELEMENT_FOLDER . "/" . $_ ; # folder with specific element of interest            
            my $tirtsd_file = $ELEMENT_FOLDER . "/" . $_ . "/" . $_ . ".tirtsd";

            # check if a manual review is already present in the README file for this element
            my $grep_res = `grep "Manual Review 1 result" $specific_element_folder/README.txt`;
            if ($grep_res) {
                warn "WARNING: A prior manual review result was found for element $_, ignoring this elment\n";
            }
            else {
                if (-e $tirtsd_file) {
                    $files{$_}[0]=$tirtsd_file;
                }
                else {
                    die "ERROR: Folder $ELEMENT_FOLDER/$_ does not have a .tirtsd file, cannot continue without this file.\n";
                }

                my $alignement_file = $ELEMENT_FOLDER . "/" . $_ . "/" . $_ . ".maf";
                if (-e $alignement_file) {
                    $files{$_}[1]=$alignement_file;
                }
                else {
                    die "ERROR: Folder $ELEMENT_FOLDER/$_ does not have a .maf file, cannot continue without this file.\n";
                }
            }
        }
    }
    unless (keys %files) {
        warn "WARNING: No files that need manual review were found\n";
    }

    my $count = 0; # used to report to the user how many elements need to be reviewed
    ## Read the .tirtsd files and process each relevant line (based on %EXAMINE_CODES)
    foreach my $element_name (keys %files) {

        $count++;
        my $review_elements = keys %files; # total number of elements to review
        print "\tElement $element_name $count of $review_elements\n";

        # Variables specific to this section
        my $TIR_b1; # left bound of TIR accepted by user
        my $TIR_b2; # right bound of TIR accpted by user
        my $TSD_size; # size of TSD accepted by user
        my $TSD_type; # if empty then it's a number othwise it's TA
        my $TIR_size; # size of TIR 

        # open the README file
        open (README, ">>$ELEMENT_FOLDER/$element_name/README.txt") or die "ERROR: Could not open or create README file $ELEMENT_FOLDER/$element_name/README.txt\n";

        # record all the relevant lines from the current .tirtsd file
        my %locs; # holds the line from from the .tirtsd file as key and priority as value
        my $filename_tirtsd = $files{$element_name}[0];
        open (INPUT, $filename_tirtsd) or die "ERROR: Cannot open file $filename_tirtsd, $!";
        <INPUT>; # skip the header
        my $i=0;
        while (my $line = <INPUT>) { # reading the individual .tirtsd file
            my @d = split " ", $line;
            $locs{$line} = $EXAMINE_CODES{$d[0]}; # store the line along with priority
            $i++;
        }
        unless ($i) { # check that at least one line was added
            die "ERROR: file $filename_tirtsd is empty, cannot continue the analysis\n";
        }

        my $menu1 = 1; # boolean, set to one until the user is done with menu 1
        while ($menu1) {

            my %alignment_sequences = fastatohash($files{$element_name}[1]); # load the existing alignment

            # Display the lines so the user can see what's available
            print "\nMENU 1\n";
            print "0) Quit this element\n\n";
            print "#) code | TIR-boundaries, # sequences with TIRs, Mean TIR length | Number of TSDs, Selected TSD\n";
            my $i=1;
            my @selections; # holds the information on the lines presented to the user

            foreach my $line (sort { $locs{$a} <=> $locs{$b} } keys %locs) {
                my @d = split " ", $line;
                my $TSD; # identity of the TSD with maximum abundance
                my $max_tsd; # highest number of TSDs on the line
            
                # determine the TSD with the highest number
                if ($d[4] > $max_tsd) { $TSD="TA"; $max_tsd =$d[4]; }
                if ($d[5] > $max_tsd) { $TSD=2; $max_tsd =$d[5];  }
                if ($d[6] > $max_tsd) { $TSD=3; $max_tsd =$d[6];  }
                if ($d[7] > $max_tsd) { $TSD=4; $max_tsd =$d[7];  }
                if ($d[8] > $max_tsd) { $TSD=5; $max_tsd =$d[8];  }
                if ($d[9] > $max_tsd) { $TSD=6; $max_tsd =$d[9];  }
                if ($d[10] > $max_tsd) { $TSD=7; $max_tsd =$d[10];  }
                if ($d[11] > $max_tsd) { $TSD=8; $max_tsd =$d[11];  }
                if ($d[12] > $max_tsd) { $TSD=9; $max_tsd =$d[12];  }
                if ($d[13] > $max_tsd) { $TSD=10; $max_tsd =$d[13];  }

                # determine the average TIR length for the current combination of sequences and locations";
                my $number_of_sequences; # total sequences in this alignment
                my $total_TIR_length; # sum of all the TIR length, used to calculate the average
                my $TIR_number; # number of sequences with TIRs
                my $average_TIR_length;

                foreach my $key (keys %alignment_sequences) {
                    my ($tir1, $tir2) = gettir($alignment_sequences{$key}, $d[1], $d[2], $MIN_TIR_SIZE, $TIR_MISMATCHES);
                    if ($tir1) {
                        $total_TIR_length += length ($tir1);
                        $TIR_number += 1;
                    }
                }

               
                if ($TIR_number) { # if TIRs have been found
                    $average_TIR_length = int($total_TIR_length / $TIR_number);
                }
                else {
                    $TIR_number = "No TIRs found";
                    $average_TIR_length = "N/A";
                }

                print "$i) $d[0] | $d[1]-$d[2], $TIR_number, $average_TIR_length | $d[4]-$d[5]-$d[6]-$d[7]-$d[8]-$d[9]-$d[10]-$d[11]-$d[12]-$d[13], $TSD\n";
                push @selections, "$d[1]\t$d[2]\t$TSD\t$average_TIR_length\n";   
                $i++;          
            }
            print "\n$i) manually enter TIRs and TSD\n"; # display menu item to enter manual coordinates
                            
            do { # read the user input until it's a number within range
                 print "Line selection: ";
                 $pkey = <STDIN>;
            } until ((looks_like_number($pkey)) and ($pkey <= $i));
            
            if ($pkey == $i) { # This means the manual selection was entered  
                my $entry_accepted; # used to know if user has finished selecting             
                do {
                    $entry_accepted = 1; # assume user will put in a correct entry, set to zero if not
                    print "Enter the coordinates manually in the form \"TIR1-left-bound TIR2-right-bound TSD-size TIR-length\" (can put \"TA\" for size, TIR-length optional)\n";
                    $pkey = <STDIN>;
                    chomp $pkey;
                    my @d = split " ", $pkey;
                    if (looks_like_number($d[0])) { $TIR_b1 = $d[0] } else { $entry_accepted = 0}
                    if (looks_like_number($d[1])) { $TIR_b2 = $d[1] } else { $entry_accepted = 0}
                    if (($d[2] eq "TA") or ($d[2] eq "ta")) {
                        $TSD_size = 2;
                        $TSD_type = "TA";
                    }
                    elsif (looks_like_number($d[2])) {
                        $TSD_size = $d[2];
                    }
                    else {
                        $entry_accepted = 0
                    }
                    if (looks_like_number($d[3])) { $TIR_size = $d[3]} else { $TIR_size = "N/A"}    
                } until ($entry_accepted);
            }
            elsif ($pkey == 0) { # the user has decided to quit
                $menu1 = 0;
            }
            else { # This means that the user selected a preset number
                my @e = split " ", $selections[$pkey-1];
                if ($e[2] eq "TA") { # adjust in case the TSD is "TA" rather than a number
                    $e[2] = 2;
                    $TSD_type = "TA";
                }
                $TIR_b1 = $e[0];
                $TIR_b2 = $e[1];
                $TSD_size = $e[2];
                $TIR_size = $e[3];
            }
            $pkey = ""; # reset the pressed key 

            if ($menu1) { # only continue if the user has not elected to quit menu 1
                # The TIRs and TSDs location have now been selected, next create an alignment to display these to the user  
                open (OUTPUT, ">", $temp_alignment_file) or die "Cannot create temporary alignment file $temp_alignment_file\n";

                foreach my $seq_name (keys %alignment_sequences) {
                    unless ($seq_name =~ /consensus/) { # avoid the line with the consensus sequence
                        print OUTPUT ">$seq_name\n";
                        ## left side sequences
                        my $left_whole_seq = substr($alignment_sequences{$seq_name}, 0, $TIR_b1);
                        $left_whole_seq =~ s/-//g; #remove gaps
                        # if there are no or few sequences, replace left TSD with space symbols
                        if ((length $left_whole_seq) < $TSD_size) {
                            $left_whole_seq = "";
                            for (my $i=0; $i<=$TSD_size; $i++) {
                                $left_whole_seq .= "s";
                            }
                        }
                        my $left_tir_seq1 = substr($left_whole_seq, -$TSD_size-1, $TSD_size);
                        my $left_tir_seq2 = substr($alignment_sequences{$seq_name}, $TIR_b1-1, $TIR_bp);

                        ## right side sequences
                        my $right_whole_seq = substr($alignment_sequences{$seq_name}, $TIR_b2, -1);
                        $right_whole_seq =~ s/-//g; #remove gaps
                        # if there are no or few sequences, replace right TSD with space symbols
                        if ((length $right_whole_seq) < $TSD_size) {
                            $right_whole_seq = "";
                            for (my $i=0; $i<=$TSD_size; $i++) {
                                $right_whole_seq .= "s";
                            }
                        }
                        my $right_tir_seq1 = substr($right_whole_seq, 0, $TSD_size);
                        my $right_tir_seq2 = substr($alignment_sequences{$seq_name}, $TIR_b2-$TIR_bp, $TIR_bp);

                        print OUTPUT $left_tir_seq1, "sss", $left_tir_seq2, "ssssssssssssssssssss", $right_tir_seq2, "sss", $right_tir_seq1, "\n";
                    }
                }

                # MENU 2 display the alignement to the user and ask for evaluation
                `aliview $temp_alignment_file`;
                if ($?) { die "Error executing: aliview $temp_alignment_file, error code $?\n"}

                my $menu2 = 1; # boolean, set to 1 until the user is done with menu 2
                my $element_rejected = 0; # boolean, set to 0 unless option "this is not an element selected", used to know which README to edit

                while ($menu2) { #keep displaying until the user ready to leave
                    print "\nMENU 2: Select what to do with this element:\n";
                    print "0) Go back to the previous menu\n";
                    if ($TSD_type eq "TA") {
                        print "1) Update the README to say this is an element with TSDs of type TA and TIRs of size $TIR_size\n";
                    }
                    else {
                        print "1) Update the README to say this as an element with TSDs of size $TSD_size and TIRs of size $TIR_size\n";
                    }
                    print "2) Update the README to say this is not an element\n";
                    print "3) Change the TSD size and/or TIR size\n";
                    print "4) Make a note in the README file\n";
                    print "5) Done reviewing this element\n";

                    do { # read the user input until it's a number within range
                        print "Line selection: ";
                        $pkey = <STDIN>;
                    } until ((looks_like_number($pkey)) and ($pkey <= 5));

                    # process the user's choice
                    if ($pkey == 0) {
                        $menu2 = 0;
                    } 
                    elsif ($pkey == 1) {
                        my $datestring = localtime(); 
                        if ($TSD_type eq "TA") {
                            print README "$datestring, Manual Review 1 result: This is an element, TSD TA, TIR size $TIR_size\n";
                        }
                        else {
                            print README "$datestring, Manual Review 1 result: This is an element, TSD $TSD_size, TIR size $TIR_size\n";
                        }
                    }
                    elsif ($pkey == 2) {
                        my $datestring = localtime(); 
                        print README "$datestring, Manual Review 1 result: This is not an element\n";
                        print REJECT "$datestring\t$element_name\tSTEP 3\tManual review of TSD and TIRs determined this is not an element\n";
                        `mv $ELEMENT_FOLDER/$element_name $ANALYSIS_FILES_OUTPUT_DIR`;
                        if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FILES_OUTPUT_DIR: error code $?\n"}
                        $element_rejected = 1;
                    }
                    elsif ($pkey == 3) {
                        my $pkey2;
                        do { # read the user input until it's a number
                            print "Enter TSD size or type \"TA\": ";
                            $pkey2 = <STDIN>;
                            chomp $pkey2;
                            if (($pkey2 eq "TA") or ($pkey2 eq "ta")) {
                                $TSD_type = "TA";
                                $pkey2 = 2;
                            }
                        } until (looks_like_number($pkey2));
                        $TSD_size = $pkey2;

                        do { # read the user input until it's a number
                            print "Enter TIR size, enter zero if size is unknown: ";
                            $pkey2 = <STDIN>;
                            chomp $pkey2;
                        } until (looks_like_number($pkey2));
                        $TIR_size = $pkey2;
                    }
                    elsif ($pkey == 4) {
                        my $datestring = localtime(); 
                        if ($element_rejected) {
                            print "Edit the file $ANALYSIS_FILES_OUTPUT_DIR/$element_name/README.txt starting with\n";
                        }
                        else {
                            print "Edit the file $ELEMENT_FOLDER/$element_name/README.txt starting with\n";
                        }
                        print "$datestring, Manual Review 1 user note: \n";
                        print "Press enter when done: ";
                        <STDIN>;
                    }
                    elsif ($pkey == 5) {
                        $menu2 = 0;
                        $menu1 = 0;
                    }
                }
            }            
        }
        close README;
    }
}
close ANALYSIS;
close REJECT;
