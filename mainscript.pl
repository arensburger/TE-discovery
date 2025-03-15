# Feb 2025 - Top-level script for TE discovery pipeline
# Inputs: List of TE sequences (FASTA) and a genome (FASTA)

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use FindBin::libs;  # this sets up that the directory lib will have all modules necessary to run the program 
use TEdiscovery;    # these are the subcripts necessary to run the pipeline
use List::UtilsBy qw(max_by);
use List::Util qw(max);

### CONSTANTS hardware
my $NUM_THREADS = 8; # number of threads to use

### CONSTANTS filtering which proteins to keep for analysis
###      using Goubert et al. as insparation and a tblastn keep hits that have the following characteristics 
my $GENOME_IDENTITY = 80; # IDENTIFYING PROTEINS, per protein, minimum percent identity between protein and genome
my $COVERAGE_RATIO = 0.5; # IDENTIFYING PROTEINS, per protein, minimum ratio of (blast match length) / (query length)
my $COPY_NUMBER = 2; # IDENTIFYING PROTEINS, minimum number of copies that hit different parts of the genome 
my $MIN_DISTANCE = 10000;   # IDENTIFYING PROTEINS, if two elements are on the same chromosome, how far they have to be, to be considered different elements
                            # NOTE: the minimum distance should bigger than the $BLAST_EXTEND variable, to avoid having the same element recorded twice        
my $BLAST_EXTEND = 2000; # IDENTIFYING PROTEINS, number of bp to extend on each side of the blast hit
my $MAX_SEQUENCE_NUMBER = 100; # ALIGNING SEQUENCES maximum number of sequences to consider, to save time
my $GAP_THRESHOLD=0.75; # REMOVING GAPS FROM ALIGNMENT, if an alignment position has this proportion or more of gaps, then remove it from the multiple sequence alignment
my $CONSLEVEL=0.60; # MAKING CONSENSUS OF SEQUENCES sequence consensus level for consensus
my $WINDOW_SIZE = 15; # MAKING CONSENSUS OF SEQUENCES size of the window looking for stretches of N's
my $SCAN_PROP = 0.5; # MAKING CONSENSUS OF SEQUENCES minimum proportion of side scan that has to be N's to be a real N full edge
my $MAX_WIN_N = 2; # MAKING CONSENSUS OF SEQUENCES maximum number of N's in the first window where the transition from N to non-N is
my $EDGE_TEST_PROPORTION = 0.05; # TESTING CONSENSUS SEQUENCES how far from the edge of the consensus do the non-gap positions have to start
my $SCAN_SIZE = 10; # IDENTIFYING TIR-TSDS number of bp to scan for TIRs from both ends
my $MIN_MATCH = 8; # IDENTIFYING TIR-TSDS minimum number of positions that match to call something a TIR
my $TIR_PROP_CUTOFF = 0.15; # IDENTIFYING TIR-TSDS proportion of elements with TIRs at which positions are reported
my $MIN_PROPORTION_SEQ_WITH_TIR=0.25; #IDENTIFYING TIRs minimum proportion of total elements for a sequence that must contain proper TIRs to be considered a candidate
my $MAX_TIR_PROPORTION=0.75; #IDENTIFYING TIRs how close to the maximum number of tirs do you have to be to qualify as a top TIR
my $MAX_END_PROPORTION=0.75; #IDENTIFYING TIRs how close to maximum proportion of sequences with identical start and stop of tir sequences you can be to a top number
my $MAX_TSD_PROPORTION=0.5; #IDENTIFYING TIRs how close to maximum number of TSDs to qualify as a top TSD sequence


### INPUTs from command line
my $INPUT_PROTEIN_SEQUENCES; # fasta formated file with input protein sequences
my $INPUT_GENOME; # fasta formated file with genome that input proteins
my $ANALYSIS_NAME; # name to give to this analysis, can be any string
my $ANALYSIS_FILES_OUTPUT_DIR; # directory where analsysis output files of the analysis are stored
my $ELEMENT_FOLDER; # directory where individual folders for each element are stored
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

## CHECK INPUTS Create output directory for analysis files if necessary
print "Preliminary steps...\n";
my $ANALYSIS_FILES_OUTPUT_DIR="./$ANALYSIS_NAME-analysis-files"; # directory to store output files of current analysis (no slash at the end), these can be destroyed when analysis is finished
if (-d $ANALYSIS_FILES_OUTPUT_DIR) {
    print "\tWARNING: Directory $ANALYSIS_FILES_OUTPUT_DIR already exists, using it to store files generated during this analysis (existing files with the same name will be overwriten)\n";
}
else {
    print "\tCreating directory $ANALYSIS_FILES_OUTPUT_DIR for storing files generated during the analysis\n";
    mkdir( $ANALYSIS_FILES_OUTPUT_DIR ) or die "Couldn't create $ANALYSIS_FILES_OUTPUT_DIR directory, $!";
}

## CHECK INPUTS Create output directory for individual elements if necessary
my $ELEMENT_FOLDER="./$ANALYSIS_NAME-elements"; # directory where the analysis of individual element will be stored
if (-d $ELEMENT_FOLDER) {
    print "\tWARNING: Directory $ELEMENT_FOLDER already exists, using it to put folders for individual elements (any current subfolders will not be overwritten and will be used for subsequent analysis)\n";
}
else {
    print "\t Creating directory $ELEMENT_FOLDER that will have subdirectories for individual elements\n";
    mkdir( $ELEMENT_FOLDER ) or die "Couldn't create $ELEMENT_FOLDER directory, $!";
}

## VARIABLES used by more than one step in the pipeline
my $blast_output_file_name = "$ANALYSIS_FILES_OUTPUT_DIR/tablastn.o"; 

## record analysis parameters
my $datestring = localtime();
open (OUTPUT,'>', "$ANALYSIS_FILES_OUTPUT_DIR/Analysis_parameters.txt") or die "ERROR: cannot create file $ANALYSIS_FILES_OUTPUT_DIR/Analysis_parameters.txt\n";
print OUTPUT "Analysis name: $ANALYSIS_NAME\n";
print OUTPUT "Date and time: $datestring\n";
print OUTPUT "Input file: $INPUT_PROTEIN_SEQUENCES\n";
print OUTPUT "Genome: $INPUT_GENOME\n\n";
print OUTPUT "Parameters, set in the script:\n";
print OUTPUT "GENOME_IDENTITY = 80; # IDENTIFYING PROTEINS, per protein, minimum percent identity between protein and genome
COVERAGE_RATIO = 0.5; # IDENTIFYING PROTEINS, per protein, minimum ratio of (blast match length) / (query length)
COPY_NUMBER = 2; # IDENTIFYING PROTEINS, minimum number of copies that hit different parts of the genome
MIN_DISTANCE = 10000;   # IDENTIFYING PROTEINS, if two elements are on the same chromosome, how far they have to be, to be considered different elements
                            # NOTE: the minimum distance should bigger than the BLAST_EXTEND variable, to avoid having the same element recorded twice        
BLAST_EXTEND = 2000; # IDENTIFYING PROTEINS, number of bp to extend on each side of the blast hit
MAX_SEQUENCE_NUMBER = 100; # ALIGNING SEQUENCES maximum number of sequences to consider, to save time
GAP_THRESHOLD=0.75; # REMOVING GAPS FROM ALIGNMENT, if an alignment position has this proportion or more of gaps, then remove it from the multiple sequence alignment
CONSLEVEL=0.60; # MAKING CONSENSUS OF SEQUENCES sequence consensus level for consensus
WINDOW_SIZE = 15; # MAKING CONSENSUS OF SEQUENCES size of the window looking for stretches of N's
SCAN_PROP = 0.5; # MAKING CONSENSUS OF SEQUENCES minimum proportion of side scan that has to be N's to be a real N full edge
MAX_WIN_N = 2; # MAKING CONSENSUS OF SEQUENCES maximum number of N's in the first window where the transition from N to non-N is
EDGE_TEST_PROPORTION = 0.05; # TESTING CONSENSUS SEQUENCES how far from the edge of the consensus do the non-gap positions have to start
SCAN_SIZE = 10; # IDENTIFYING TIR-TSDS number of bp to scan for TIRs from both ends
MIN_MATCH = 8; # IDENTIFYING TIR-TSDS minimum number of positions that match to call something a TIR
TIR_PROP_CUTOFF = 0.15; # IDENTIFYING TIR-TSDS proportion of elements with TIRs at which positions are reported
MIN_PROPORTION_SEQ_WITH_TIR=0.25; #IDENTIFYING TIRs minimum proportion of total elements for a sequence that must contain proper TIRs to be considered a candidate
MAX_TIR_PROPORTION=0.75; #IDENTIFYING TIRs how close to the maximum number of tirs do you have to be to qualify as a top TIR
MAX_END_PROPORTION=0.75; #IDENTIFYING TIRs how close to maximum proportion of sequences with identical start and stop of tir sequences you can be to a top number
MAX_TSD_PROPORTION=0.5; #IDENTIFYING TIRs how close to maximum number of TSDs to qualify as a top TSD sequence\n\n";
print OUTPUT "TIR-TSD test for a specific TIR location:\n";
print OUTPUT "These tests are coded by 0 and 1's in the file ending in .tirtsd (if TIRs are found for this element)\n";
print OUTPUT "Test1: Do sequences with this TIRs make up at least $MIN_PROPORTION_SEQ_WITH_TIR of all the sequences for this element?\n";
print OUTPUT "Test2: Is this one of the most abundant TIRs (i.e. it's at least $MAX_TIR_PROPORTION times the maximum number of TIRs seen for this element)?\n";
print OUTPUT "Test3: TIR copies should generally start with the same nucleotide, and end with the same nucleotide. Is number of sequences for this TIR at least $MAX_END_PROPORTION of the highest observed TIRs that start and end with the same nucleotides?\n";
print OUTPUT "Test4: TIRs should have many TSDs associated with it, does this TIR have at least $MAX_TSD_PROPORTION of the maximum of number of TSDs seen for this element?\n";
close OUTPUT;

### PIPELINE STEP 1 identify proteins that match the genome with parameters specified above under "CONSTANTS"
###     The output is a list of proteins for further analysis recorded in the file $output_file_name
my $step_number = 1;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## VARIABLES, variable for this step
    my %protein_ids; # holds the id the input proteins that passed the filtering tests as key and the number of copies that passed the test as values

    ## Excute the tblastn search
   `tblastn -query $INPUT_PROTEIN_SEQUENCES -db $INPUT_GENOME -outfmt "6 qseqid sseqid sstart send pident length qlen" -out $blast_output_file_name -num_threads $NUM_THREADS`;
   if ($?) { die "ERROR executing tblastn, stopping analysis (hint: was the genome formated with makeblastdb?): error code $?\n"}

    ## Inspired by the Goubert et al. protocol, filter elements that 1) have >= 80% identity to genome, 2) have 50% length of the query, 3) are found at multiple locations
    my %candidate_protein; # hash with protein name as key and string with chromosome and middle location of element on that chromsome

    open (INPUT, "$blast_output_file_name") or die "ERROR: Cannot open file $blast_output_file_name\n";
    while (my $line = <INPUT>) {
        my $gi=0; # boolean, set to zero until the genome identity test is passed
        my $cr=0; # boolean, set to zero until the coverage ratio test is passed
        my $md=1; # boolean, set to one unless minimum distance test is failed

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

    # Filter out elements that have too few copy numbers
    # (making this a separate step so the code is more modular, rather then incoroporating it into the next step)
    # record any elements that were discarded at this point
    open (OUTPUT, '>', "$ANALYSIS_FILES_OUTPUT_DIR/Input_sequences_too_few_blast_hits.txt") or die "ERROR, cannot create output file $ANALYSIS_FILES_OUTPUT_DIR/Input_sequences_too_few_blast_hits.txt\n";
    print OUTPUT "These input sequences were not analyzed because the BLAST analysis returned too few hits to the genome\n";
    print OUTPUT "Input_name\tMinimum_BLAST_hits\tObserved_BLAST_hits\n";
    foreach my $prot_name (keys %protein_ids) {
        unless ($protein_ids{$prot_name} >= $COPY_NUMBER) {
            print OUTPUT "$prot_name\t$COPY_NUMBER\t$protein_ids{$prot_name}\n";
            delete $protein_ids{$prot_name};
        }
    }
    close OUTPUT;

    # Create individual directories for each element
    my $i=0; # counts the number of output lines, to check if it's zero
    foreach my $prot_name (keys %protein_ids) {
        my $folder_name = "$ELEMENT_FOLDER/$prot_name";
        if (-d $folder_name) {
            print "\tElement folder $folder_name already exists, using it without modification to existing files\n";
        }
        else {
            mkdir( $folder_name ) or die "Couldn't create $folder_name directory, $!";
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
        open (INPUT, $blast_output_file_name) or die "ERROR: cannot open file blast file output $blast_output_file_name\n";
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
     
       `mafft --quiet --thread -$NUM_THREADS $extended_fasta_name > $aligned_sequences`;
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
            }
        }

        # output the alignment file with positions removed, called .maf
        open (OUTPUT, '>', "$ELEMENT_FOLDER/$element_name/$element_name.maf") or die "ERROR: Cannot create file $element_name.maf, $!\n";
        foreach my $name (keys %seqrmg) {
            print OUTPUT ">$name\n$seqrmg{$name}\n";
            $ali_trimmed_length = length $seqrmg{$name}; # legnth of the alignment after trimming
        }
        my $datestring = localtime(); # update the README file of new file created
        print README "$datestring, File $element_name.maf is the alignment of sequences with positions containing $GAP_THRESHOLD proportion of gaps removed\n";
        close OUTPUT;

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

        # if both a left and right transition are found, then create the trimmed consensus sequence
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

        # Test to see if this consensus sequence passes the "edge test". The trimmed consensus should not be too close 
        # to the edge of the alignment (this would suggest that the whole region, not just an element inside it, is conserved)
        my $edge_test; # boolean, will be set to 1 if the consensus is not too close to the edge
        if (($ltrans > ($EDGE_TEST_PROPORTION * (length $trimmed_conseq))) and ($rtrans > ((1 - $EDGE_TEST_PROPORTION) * (length $trimmed_conseq)))) {
            $edge_test = 1;
        }

        # update the README file
        my $datestring = localtime(); 
        if ($edge_test and $trimmed_conseq) {
            print README "$datestring, a trimmed consensus sequence was created that passed the \"edge test\"\n";
        }
        elsif ($trimmed_conseq) { # get here if a consenus was found but it failed the edge test
            print README "$datestring, a trimmed consensus sequence was created but it failed the \"edge test\", stopping the analysis here\n";
        }
        else { # get here if a no abrupt changes in the consensus sequences was detected
            print README "$datestring, the consensus sequence showed no transitions into an element, stopping the analysis here\n";
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
		            my $number_of_tirs_found; # nubmer of sequences that match the TIR criteria
                    my %tir_first_and_last_bases; # first and last set of bases of tir as key and abundance as value
                    my %tsds_found; # keys is TSD type "TA", "2", ... "10" and key is number of TSDs found
                    foreach my $sequence_name (keys %seqrmg) {
                        my ($tir_found, $first_base, $last_base) = gettir($seqrmg{$sequence_name}, $i, $j); # figure if this sequences has a tir at these positions and if so, report first and last nucleotide
                        $number_of_tirs_found += $tir_found;
                        if ($tir_found) {
                            my $bases = $first_base . $last_base;
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
                    my $most_abundant_tir_seq;
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

            # report if no TIR-TSD combnations have been found
            unless (scalar @tsd_tir_combinations) {
                my $datestring = localtime(); 
                print README "$datestring, no acceptable TIR-TSD combinations have been found, stopping the analysis here.\n";
            }

            # go through the element and identify those candidate locations that pass the tests for TIR-TSD combinatations
            my @candidates_tsdtir; # locations and tsd numbers of candidate locations the tests along with success code 
            my $i; # counter of the number of locations that passed all the tests
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

                # record the results for this candidate location
                my $success_code = $min_prop_seq_wtir . $top_tir_number . $top_end_proportions . $top_number_tsds;
                push @candidates_tsdtir, "$success_code\t$d[0]\t$d[1]\t$d[2]\t$d[4]\t$d[5]\t$d[6]\t$d[7]\t$d[8]\t$d[9]\t$d[10]\t$d[11]\t$d[12]\t$d[13]";
                if ($success_code eq "1111") {
                    $i++;
                }
            }

            # Record the candidate locations and positions
            my $number_of_candidates = scalar @candidates_tsdtir;
            my $datestring = localtime(); 
            print README "$datestring, Examined $number_of_candidates candidate locations for TSD-TIR combinations, found $i locations that passed all the tests (see Analysis_paramters.txt file for details of the tests)\n";
            open (OUTPUT,'>', "$ELEMENT_FOLDER/$element_name/$element_name.tirtsd") or die "Error: cannot create file $ELEMENT_FOLDER/$element_name/$element_name.tirtsd, $!\n";
            print OUTPUT "# success_code\tloc1\tloc2\tnumber_of_tirs\tnumber_tsds_TA_through_10\n";
            foreach my $c (@candidates_tsdtir) {
                print OUTPUT "$c\n";
            }
        }     
        close README;
    }    
}

sub gettir {
	my ($seq, $loc1, $loc2) = @_;

	### load the sequence into memory
	my $sequence = substr($seq,$loc1-1,$loc2-$loc1+1); # DNA sequence of the consensus
	### get the ends into string Variables
	my $s1 = substr ($sequence, 0, $SCAN_SIZE);
	my $s2 = substr ($sequence, -$SCAN_SIZE, $SCAN_SIZE);
	# reverse complement seq2
	my $s2rc = rc($s2);
	### count the matches
	my $matches;
	for (my $i=0; $i < length $s1; $i++){
		my $char_s1 = substr ($s1, $i, 1);
		my $char_s2 = substr ($s2rc, $i, 1);

		my $countN_char1 = () = $char_s1 =~ /N|n|-/; # count forbiden characters
		my $countN_char2 = () = $char_s2 =~ /N|n|-/;
		if (($i==0) and ($countN_char1 or $countN_char2)) { # check if the TIR starts with forbidden character
			return (0);
		}

		unless ($countN_char1 or $countN_char2) {
			if ($char_s1 eq $char_s2) {
				$matches++;
			}
		}
	}
	### evaluate the results
	if ($matches >= $MIN_MATCH) {
        my $c1 = substr($s1, 0, 3);
        my $c2 = substr($s2, -3, 3);
		return (1, $c1, $c2); # tir found, report that it was found and bases
	}
	else {
		return (0, "", ""); # tir not found
	}
}

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
    
    if (($type eq "TA") and ($c1 eq "ta") and ($c2 eq "ta")) {
        return (1); # found a TA TSD
    }
    elsif (($type eq "2") and ($c1 and $c2) and ($c1 eq $c2)) { # check that both are equal and not 0, 0 means the sequence contained an "n" charcater
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

#reverse complement
sub rc {
    my ($sequence) = @_;
    $sequence = reverse $sequence;
    $sequence =~ tr/ACGTRYMKSWacgtrymksw/TGCAYRKMWStgcayrkmws/;
    return ($sequence);
}


#takes a string, converts it to lower case and checks if it contain an "n"
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
