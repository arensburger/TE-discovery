# Feb 2025 - Top-level script for TE discovery pipeline
# Inputs: List of TE sequences (FASTA) and a genome (FASTA)

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use FindBin::libs;  # this sets up that the directory lib will have all modules necessary to run the program 
use TEdiscovery;    # these are the subcripts necessary to run the pipeline
use List::UtilsBy qw(max_by);
use List::Util qw(max);

### CONSTANTS for all steps
my $NUM_THREADS = 8; # number of threads to use

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
    print "\tWARNING: Directory $ANALYSIS_FILES_OUTPUT_DIR already exists, using it to store files generated during this analysis (existing files with the same name will be overwriten, but will not overwrite folders, that will crash the script)\n";
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
    print "\tCreating directory $ELEMENT_FOLDER that will have subdirectories for individual elements\n";
    mkdir( $ELEMENT_FOLDER ) or die "Couldn't create $ELEMENT_FOLDER directory, $!";
}

## VARIABLES used by more than one step in the pipeline
my $blast_output_file_name = "$ANALYSIS_FILES_OUTPUT_DIR/tblastn.o"; 

## Create and record start parameters in file
my $datestring = localtime();
my $analysis_parameters_file_name = "$ANALYSIS_FILES_OUTPUT_DIR/Analysis_parameters.txt";
if (-f $analysis_parameters_file_name) {
    print "\tWARNING: Analysis parameters file $analysis_parameters_file_name already exists, appending a new set of parameters to that file\n";
    open (ANALYSIS,'>>', $analysis_parameters_file_name) or die "ERROR: cannot open file $analysis_parameters_file_name\n";
}
else {
    open (ANALYSIS,'>', $analysis_parameters_file_name) or die "ERROR: cannot open file $analysis_parameters_file_name\n";
}
print ANALYSIS "Analysis name: $ANALYSIS_NAME\n";
print ANALYSIS "Date and time: $datestring\n";
print ANALYSIS "Input file: $INPUT_PROTEIN_SEQUENCES\n";
print ANALYSIS "Genome: $INPUT_GENOME\n\n";
print ANALYSIS "Parameters, set in the script:\n";

## Create file to record rejected sequences
my $rejection_file_name = "$ANALYSIS_FILES_OUTPUT_DIR/Rejected_sequences.txt";
if (-f $rejection_file_name) {
    print "\tWARNING: Rejected sequences text already exists, adding new entries to this file\n";
    open (REJECT, '>>', $rejection_file_name) or die "ERROR, cannot open output file $rejection_file_name\n";
}
else {
    open (REJECT, '>', $rejection_file_name) or die "ERROR, cannot create output file $rejection_file_name\n";
}

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
print ANALYSIS "STEP2: MAX_TSD_PROPORTION = $MAX_TSD_PROPORTION\n";
print ANALYSIS "STEP2: Success code, first number: Is the proportion of sequences with TIRs higher than MIN_PROPORTION_SEQ_WITH_TIR?\n";
print ANALYSIS "STEP2: Success code, second number: Does this candidate have one of the highest number of TIRs?\n";
print ANALYSIS "STEP2: Success code, third number: Do the TIRs tend to start and end with the same bases?\n";
print ANALYSIS "STEP2: Success code, first number: Does this candidate have one of the highest number of TSDs?\n";

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

            # # report if no TIR-TSD combnations have been found
            # unless (scalar @tsd_tir_combinations) {
            #     my $datestring = localtime(); 
            #     print README "$datestring, no acceptable TIR-TSD combinations have been found, stopping the analysis here.\n";
            #     print REJECT "$datestring\t$element_name\tSTEP 2\tno TIR-TSD found\n";
            #     `mv $ELEMENT_FOLDER/$element_name $ANALYSIS_FILES_OUTPUT_DIR`;
            #     if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FILES_OUTPUT_DIR: error code $?\n"}
            # }

            # go through the element and identify those candidate locations that pass the tests for TIR-TSD combinatations
            my @candidates_tsdtir; # locations and tsd numbers of candidate locations the tests along with success code 
            my $i=0; # counter of the number of locations that passed all the tests
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
                #print "4 is $d[4], 5 is $d[5], 6 is $d[6], 7 is $d[7], 8 is $d[8], 9 is $d[9], 10 is $d[10], 11 is $d[11], 12 is $d[12], 13 is $d[13]\n";
                push @candidates_tsdtir, "$success_code\t$d[0]\t$d[1]\t$d[2]\t$d[4]\t$d[5]\t$d[6]\t$d[7]\t$d[8]\t$d[9]\t$d[10]\t$d[11]\t$d[12]\t$d[13]";
                if ($success_code eq "1111") {
                    $i++;
                }
            }

            # Record the candidate locations and positions
            my $number_of_candidates = scalar @candidates_tsdtir;
            my $datestring = localtime(); 
            print README "$datestring, Examined $number_of_candidates candidate locations for TSDs and TIRs, found $i locations that passed all the tests\n";
            open (OUTPUT,'>', "$ELEMENT_FOLDER/$element_name/$element_name.tirtsd") or die "Error: cannot create file $ELEMENT_FOLDER/$element_name/$element_name.tirtsd, $!\n";
            print OUTPUT "# success_code\tloc1\tloc2\tnumber_of_tirs\tnumber_tsds_TA_through_10\n";
            foreach my $c (@candidates_tsdtir) {
                print OUTPUT "$c\n";
            }
            my $datestring = localtime(); 
            print README "$datestring, File $element_name.tirtsd has the number of all TSDs and TIRs along with number of TSDs and codes for tests\n";
        }     
        close README;
    }    
}

### PIPELINE STEP 3 
### Present the elements to the user for manual review
### CONSTANTS applicable only for STEP 23
my %EXAMINE_CODES=("1111" => 1, "1101" => 2); # success codes to examine as key and priority as value
my $TIR_bp = 30; # how many bp to display on the TIR side

my $step_number = 3;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    my $pkey; # pressed key, used for input from user
    ## read and check the inputs
    my $ELEMENT_FOLDER="./$ANALYSIS_NAME-elements"; # don't need this when merge
    unless ($ANALYSIS_NAME) {
        die "ERROR: STEP $step_number requires an analysis name\n.";
    }

    ## make a list of files to analyze, put them into %files
    my %files; # holds the element folder path as key and path to the .tirtsd and the .maf files as array of values
    opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
    while (readdir $dh) {

         ## need to check the README file here to make sure this element has not already been reviewed manually

        unless ($_ =~ /^\./) { # prevents reading invisible files or the . and .. files
            my $tirtsd_file = $ELEMENT_FOLDER . "/" . $_ . "/" . $_ . ".tirtsd";
            if (-e $tirtsd_file) {
                $files{$ELEMENT_FOLDER}[0]=$tirtsd_file;
            }
            else {
                die "ERROR: Folder $ELEMENT_FOLDER/$_ does not have a .tirtsd file, this should not happen.\n";
            }

            my $alignement_file = $ELEMENT_FOLDER . "/" . $_ . "/" . $_ . ".maf";
            if (-e $alignement_file) {
                $files{$ELEMENT_FOLDER}[1]=$alignement_file;
            }
            else {
                die "ERROR: Folder $ELEMENT_FOLDER/$_ does not have a .maf file, this should not happen.\n";
            }
        }
    }
    unless (keys %files) {
        die "ERROR: No files that need manual review were found\n";
    }

    ## Read the .tirtsd files and process each relevant line (based on %EXAMINE_CODES)
    foreach my $element_path (keys %files) {

        # Variables specific to this section
        my $TIR_b1; # left bound of TIR accepted by user
        my $TIR_b2; # right bound of TIR accpted by user
        my $TSD_size; # size of TSD accepted by user

        # record all the relevant lines from the current .tirtsd file
        my %locs; # holds the line from from the .tirtsd file as key and priority as value
        my $filename_tirtsd = $files{$element_path}[0];
        open (INPUT, $filename_tirtsd) or die "ERROR: Cannot open file $filename_tirtsd, $!";
        while (my $line = <INPUT>) { # reading the individual .tirtsd file
            my @d = split " ", $line;
            if (exists $EXAMINE_CODES{$d[0]}) {
                $locs{$line} = $EXAMINE_CODES{$d[0]}; # store the line along with priority
            }
        }

        my $process_element = 1; # boolean, set to one until this element has been processed
        while ($process_element) {

            my %alignment_sequences = fastatohash($files{$ELEMENT_FOLDER}[1]); # load the existing alignment

            # Count the number of lines to display
            my $number_of_lines = keys %locs;
            
            if ($number_of_lines == 0) { # leave if there's nothing to do
                print "No TIR-TSD combination lines were found in file $filename_tirtsd\n";
                $process_element = 1;
            } 
            else { # there are elements to process, get the information out of the user

                # Display the lines so the user can see what's available
                print "Lines from $filename_tirtsd\n";
                print "code, TIR-boundaries, TSD numbers\n";
                foreach my $line (sort { $locs{$a} <=> $locs{$b} } keys %locs) { # sort so high priority ones are first
                    my @d = split " ", $line;
                    print "$d[0], $d[1]-$d[2], $d[4]-$d[5]-$d[6]-$d[7]-$d[8]-$d[9]-$d[10]-$d[11]-$d[12]-$d[13]\n";
                }
                print "\n";

                # Ask the user to pick a combination of TIR boundaries and TSD
                print "Select TIRs and TSDs below, or manually enter a location\n";
                print "Selection) code, TIR-boundaries, TSD type, Number of TIRS, Average TIR length\n";
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
                        print "$i) $d[0], $d[1]-$d[2], $TSD, $TIR_number, $average_TIR_length\n";  
                    }
                    else {
                        print "$i) $d[0], $d[1]-$d[2], $TSD, no TIRs found\n";
                    }
                    push @selections, "$d[1]\t$d[2]\t$TSD\n";   
                    $i++;          
                }
                print "$i) manually enter TIRs and TSD\n";
                               
                while (($pkey == 0) or ($pkey > $i)) { # read the user input
                    print "Line selection [1]: ";
                    $pkey = <STDIN>;
                }

                if ($pkey == $i) { # This means the manual selection was entered
                    my $entry_rejected = 1; # set to one until the manual entry has been accepted
                    while ($entry_rejected) {
                        print "Enter the coordinates manually in the form \"TIR1 TIR2 TSD_size\"\n";
                        $pkey = <STDIN>;
                        if ($pkey =~ /(\d+)\s(\d+)\s(\d+)/) {
                            $TIR_b1 = $1;
                            $TIR_b2 = $2;
                            $TSD_size = $3;
                            $entry_rejected = 0;
                        }
                        else {
                            print "ERROR: Did not recognize this manual entry\n";
                        }
                    }
                }
                else { # This means that the user selected a preset number
                    my @e = split " ", $selections[$pkey-1];
                    $TIR_b1 = $e[0];
                    $TIR_b2 = $e[1];
                    $TSD_size = $e[2];
                }
                $pkey = ""; # reset the pressed key 
            } 

            # The TIRs and TSDs locations have been selected, next create an alignment to display these to the user  
 #           my %alignment_sequences = fastatohash($files{$ELEMENT_FOLDER}[1]); # load the existing alignment
            my $display_alignment_file = File::Temp->new(UNLINK => 1, SUFFIX => '.fa' ); # file with alignment that will be diplayed to the user
            foreach my $seq_name (keys %alignment_sequences) {
                unless ($seq_name =~ />consensus/) { # avoid the line with the consensus sequence
                    print $display_alignment_file ">$seq_name\n";
                    ## left side sequences
                    my $left_whole_seq = substr($alignment_sequences{$seq_name}, 0, $TIR_b1);
                    $left_whole_seq =~ s/-//g; #remove gaps
                    # if there are no or few sequences, replace left TSD with space symbols
                    if ((length $left_whole_seq) < $TSD_size) {
                        $left_whole_seq = "";
                        for (my $i=0; $i<$TSD_size; $i++) {
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
                        for (my $i=0; $i<$TSD_size; $i++) {
                            $right_whole_seq .= "s";
                        }
                    }
                    my $right_tir_seq1 = substr($right_whole_seq, 0, $TSD_size);
                    my $right_tir_seq2 = substr($alignment_sequences{$seq_name}, $TIR_b2-$TIR_bp, $TIR_bp);

                    print $display_alignment_file $left_tir_seq1, "sss", $left_tir_seq2, "ssssssssssssssssssss", $right_tir_seq2, "sss", $right_tir_seq1, "\n";
                }
            }
`cp $display_alignment_file /home/peter/Desktop/ali.fa`;            
exit;

        }
    }
}

close ANALYSIS;
close REJECT;
