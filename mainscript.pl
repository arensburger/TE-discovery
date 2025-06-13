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
use Term::ANSIColor;

### INPUTs from command line, top level variables are in uppercase
my $INPUT_PROTEIN_SEQUENCES; # fasta formated file with input protein sequences
my $TBLASTN_FILE; # name of the file containing the out put of the tblastn
my $INPUT_GENOME; # fasta formated file with genome that input proteins
my $ANALYSIS_FOLDER; # name of folder to store analysis files into
my $ELEMENT_FOLDER; # directory where individual folders for each element are stored
my $REJECTED_ELEMENTS_FOLDER = "Rejected_elements"; # name of folder that contains files for elements have been reviewed and rejected
my $START_STEP; # analysis step to start at
my $END_STEP; # analysis step to end at,
my $SHOW_HELP; # call for help 

### CHECK INPUTS Read and check that the inputs have been provided
GetOptions(
	'p:s'   => \$INPUT_PROTEIN_SEQUENCES,
	't:s'   => \$TBLASTN_FILE,
	'g:s'   => \$INPUT_GENOME,
    'n:s'   => \$ANALYSIS_FOLDER,
    'e:s'   => \$ELEMENT_FOLDER,
    'a:i'   => \$START_STEP,
    'b:i'   => \$END_STEP,
    'h'     => \$SHOW_HELP,
);

## CHECK INPUTS if help was called (this could probably be improved)
if ($SHOW_HELP) {
    print "Description and input help can be found at https://github.com/arensburger/TE-discovery\n";
    exit;
}

## CHECK INPUTS Validate required input files
unless ($ANALYSIS_FOLDER and $ELEMENT_FOLDER) {
    die "usage: perl mainscript.pl <-n folder name (full path) to store analysis files REQUIRED> <-e folder name to store elements REQUIRED> <-h for more help>\n";
}

## CHECK INPUTS set the start and end steps
if ($START_STEP) {
    unless ($END_STEP) {
        $END_STEP = $START_STEP;
    }
}
else { 
    $START_STEP = 1;
    $END_STEP = 100000; # set to a very high number;
}

## CHECK INPUTS Create output directory for analysis files if necessary
#$ANALYSIS_FOLDER = fixdirname($ANALYSIS_FOLDER);
if (-d $ANALYSIS_FOLDER) {
    warn "Directory $ANALYSIS_FOLDER already exists\n";
}
else {
    print "Creating directory $ANALYSIS_FOLDER for storing files generated during the analysis\n";
    `mkdir $ANALYSIS_FOLDER`;
    if ($?) { die "ERROR creating directory: error code $?\n"}

}

## CHECK INPUTS Create output directory for elements that have been rejected
my $reject_folder_path = $ANALYSIS_FOLDER . "/" . $REJECTED_ELEMENTS_FOLDER;
if (-d $reject_folder_path) {
    warn "Directory $reject_folder_path already exists\n";
}
else {
    print "Creating directory $reject_folder_path for storing rejected elements\n";
    `mkdir $reject_folder_path`;
    if ($?) { die "ERROR creating directory: error code $?\n"}
}

## CHECK INPUTS Create output directory for individual elements if necessary
#$ELEMENT_FOLDER = fixdirname($ELEMENT_FOLDER);
if (-d $ELEMENT_FOLDER) {
    warn "Directory $ELEMENT_FOLDER already exists\n";
}
else {
    print "Creating directory $ELEMENT_FOLDER that will have subdirectories for individual elements\n";
    `mkdir $ELEMENT_FOLDER`;
    if ($?) { die "ERROR creating directory: error code $?\n"}
}

## Create or open files to store analysis parameters, and to store rejected sequences. 
my $analysis_parameters_file_name = "$ANALYSIS_FOLDER/Analysis_parameters.txt"; # file to record parameters
my $rejection_file_name = "$ANALYSIS_FOLDER/Rejected_sequences.txt"; # file to store rejected sequences, and why
open (ANALYSIS,'>>', $analysis_parameters_file_name) or die "ERROR: cannot open file $analysis_parameters_file_name\n"; # create or append to file
open (REJECT, '>>', $rejection_file_name) or die "ERROR, cannot create output file $rejection_file_name\n"; # create or append to file

my $datestring = localtime();
print ANALYSIS "\nSTARTING ANALYSIS on $datestring\n";

my $BLAST_OUTPUT_FILE_NAME = "$ANALYSIS_FOLDER/tblastn.o"; # default name and location unless a file is provided

### PIPELINE STEP 1 identify proteins that match the genome with parameters specified above under
###     The output is a list of proteins for further analysis recorded in the file $output_file_name
### CONSTANTS applicable to this step only (also record these in the file)
my $GENOME_IDENTITY = 80; # IDENTIFYING PROTEINS, per protein, minimum percent identity between protein and genome
my $COVERAGE_RATIO = 0.5; # IDENTIFYING PROTEINS, per protein, minimum ratio of (blast match length) / (query length)
my $COPY_NUMBER = 2; # IDENTIFYING PROTEINS, minimum number of copies that hit different parts of the genome 
my $MIN_DISTANCE = 10000;   # IDENTIFYING PROTEINS, if two elements are on the same chromosome, how far they have to be, to be considered different elements
                            # NOTE: the minimum distance should bigger than the $BLAST_EXTEND variable, to avoid having the same element recorded twice    
my $NUM_THREADS = `nproc --all`;# determine the number of processors on the current machine
if ($?) { warn "WARNING could not determine the number of cores automatically, defaulting to 8\n"; $NUM_THREADS=8}
chomp $NUM_THREADS; 

my $step_number = 1;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

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
            die ("ERROR, running STEP 1 requires either than both -p and -g parameters are set or that -t is set\n");
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
            print REJECT "$datestring $dupseq overlaps with $topseq and was taken out of the analysis at STEP $step_number \n";
        }

        # run tblastn
        `tblastn -query $protein_file_no_redudants -db $INPUT_GENOME -outfmt "6 qseqid sseqid sstart send pident length qlen" -out $BLAST_OUTPUT_FILE_NAME -num_threads $NUM_THREADS`;
        if ($?) { die "ERROR executing tblastn, stopping analysis (hint: was the genome formated with makeblastdb?): error code $?\n"}
        print ANALYSIS "\ttblastn was executed, the output is in file $BLAST_OUTPUT_FILE_NAME\n";
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

my $MAX_SEQUENCE_NUMBER = 100; # ALIGNING SEQUENCES maximum number of sequences to consider, to save time
my $GAP_THRESHOLD=0.75; # REMOVING GAPS FROM ALIGNMENT, if an alignment position has this proportion or more of gaps, then remove it from the multiple sequence alignment
my $CONSLEVEL=0.60; # MAKING CONSENSUS OF SEQUENCES sequence consensus level for consensus
my $WINDOW_SIZE = 15; # MAKING CONSENSUS OF SEQUENCES size of the window looking for stretches of N's
my $SCAN_PROP = 0.5; # MAKING CONSENSUS OF SEQUENCES minimum proportion of side scan that has to be N's to be a real N full edge
my $MAX_WIN_N = 2; # MAKING CONSENSUS OF SEQUENCES maximum number of N's in the first window where the transition from N to non-N is
my $EDGE_TEST_PROPORTION = 0.05; # TESTING CONSENSUS SEQUENCES how far from the edge of the consensus do the non-gap positions have to start
my $MIN_TIR_SIZE = 10; # IDENTIFYING TIR-TSDS smallest allowable size for the TIR
my $TIR_MISMATCHES = 2; # IDENTIFYING TIR-TSDS maximum number of mismatches allowed between two TIRs
my $TIR_PROP_CUTOFF = 0.15; # IDENTIFYING TIR-TSDS proportion of elements with TIRs at which positions are reported
my $MIN_PROPORTION_SEQ_WITH_TIR=0.25; #IDENTIFYING TIRs minimum proportion of total elements for a sequence that must contain proper TIRs to be considered a candidate
my $MAX_TIR_PROPORTION=0.75; #IDENTIFYING TIRs how close to the maximum number of tirs do you have to be to qualify as a top TIR
my $MAX_END_PROPORTION=0.75; #IDENTIFYING TIRs how close to maximum proportion of sequences with identical start and stop of tir sequences you can be to a top number
my $MAX_TSD_PROPORTION=0.5; #IDENTIFYING TIRs how close to maximum number of TSDs to qualify as a top TSD sequence
my %EXAMINE_CODES=("1111" => 1, "1101" => 2); # success codes to examine as key and priority as value

my $step_number = 2;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## update the analysis file with what is being done and paramter values
    print ANALYSIS "Running STEP 2\n";
    print ANALYSIS "\tBLAST_EXTEND = $BLAST_EXTEND\n";
    print ANALYSIS "\tMAX_SEQUENCE_NUMBER = $MAX_SEQUENCE_NUMBER\n";
    print ANALYSIS "\tGAP_THRESHOLD = $GAP_THRESHOLD\n";
    print ANALYSIS "\tCONSLEVEL = $CONSLEVEL\n";
    print ANALYSIS "\tWINDOW_SIZE = $WINDOW_SIZE\n";
    print ANALYSIS "\tSCAN_PROP = $SCAN_PROP\n";
    print ANALYSIS "\tMAX_WIN_N = $MAX_WIN_N\n";
    print ANALYSIS "\tEDGE_TEST_PROPORTION = $EDGE_TEST_PROPORTION\n";
    print ANALYSIS "\tMIN_TIR_SIZE = $MIN_TIR_SIZE\n";
    print ANALYSIS "\tTIR_MISMATCHES = $TIR_MISMATCHES\n";
    print ANALYSIS "\tTIR_PROP_CUTOFF = $TIR_PROP_CUTOFF\n";
    print ANALYSIS "\tMIN_PROPORTION_SEQ_WITH_TIR = $MIN_PROPORTION_SEQ_WITH_TIR\n";
    print ANALYSIS "\tMAX_TIR_PROPORTION = $MAX_TIR_PROPORTION\n";
    print ANALYSIS "\tMAX_END_PROPORTION = $MAX_END_PROPORTION\n";
    print ANALYSIS "\t%EXAMINE_CODES=(\"1111\" => 1, \"1101\" => 2)\n";
    print ANALYSIS "\tMAX_TSD_PROPORTION = $MAX_TSD_PROPORTION\n";

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
            `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
            if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $reject_folder_path: error code $?\n"}
        }
        else { # get here if no abrupt changes in the consensus sequences was detected
            print README "$datestring, the consensus sequence showed no transitions into an element, stopping the analysis here\n";
            print REJECT "$datestring\t$element_name\tSTEP 2\tNo transitions to element observed in alignment file\n";
            `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
            if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $reject_folder_path: error code $?\n"}
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
                `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
                if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FOLDER: error code $?\n"}
            }
        }     
        close README;
    }    
}

### PIPELINE STEP 3 
### Present the elements to the user for manual review
### CONSTANTS applicable only for STEP 3
my $TIR_bp = 30; # how many bp to display on the TIR side

my $step_number = 3;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

     ## update the analysis file with what is going on
    print ANALYSIS "Running STEP 3\n";
    print ANALYSIS "\tTIR_bp = $TIR_bp\n";

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
                print "\tElement $_ has already been manually reviewed, ignoring\n";
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
        print color('bold blue');
        print "\nElement $element_name $count of $review_elements\n";
        print color('reset');

        # Variables specific to this section
        my $TIR_b1; # left bound of TIR accepted by user
        my $TIR_b2; # right bound of TIR accpted by user
        my $TSD_size; # size of TSD accepted by user
        my $TSD_type; # if empty then it's a number othwise it's TA
        my $TIR_size; # size of TIR 

        # check the README file for any previous manual review notes and display them
        open (README, "$ELEMENT_FOLDER/$element_name/README.txt") or die "ERROR: Could not open or create README file $ELEMENT_FOLDER/$element_name/README.txt\n";
        my $prior_notes; 
        while (my $line = <README>) {
            if ($line =~ /Manual\sReview\s1\suser\snote/) {
                $prior_notes .= $line;
            }
        } 
        if ($prior_notes) {
            print "\nPRIOR REVIEW NOTES FOR ELEMENT $element_name\n$prior_notes";
        }
        close README;

        # open the README file for writing
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
            print color('magenta');
            print "\nMENU 1\n";
            print color('reset');
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
            print "$i) manually enter TIRs and TSD\n"; # display menu item to enter manual coordinates
                            
            do { # read the user input until it's a number within range
                 print "Line selection: ";
                 $pkey = <STDIN>;
            } until ((looks_like_number($pkey)) and ($pkey <= $i));
            
            if ($pkey == $i) { # This means the manual selection was entered  
                my $entry_accepted; # used to know if user has finished selecting             
                do {
                    $entry_accepted = 0; # assume user will put in a correct entry, set to zero if not
                    print "Enter left alignement coordinate: ";
                    $pkey = <STDIN>; chomp $pkey;
                    $TIR_b1 = $pkey;
                    print "Enter right alignement coordinate: ";
                    $pkey = <STDIN>; chomp $pkey;
                    $TIR_b2 = $pkey;
                    print "Enter TSD size (can enter TA): ";
                    $pkey = <STDIN>; chomp $pkey;
                    $TSD_size = $pkey;
                    print "Enter TIR size (can enter 0 if you don't want to enter a size): ";
                    $pkey = <STDIN>; chomp $pkey;
                    $TIR_size = $pkey;
                    if (looks_like_number($TIR_b1) and looks_like_number($TIR_b2) and looks_like_number($TIR_size)) {
                        if (($TSD_size eq "TA") or ($TSD_size eq "ta")) {
                            $TSD_size = 2;
                            $TSD_type = "TA";
                            $entry_accepted=1;
                        }
                        elsif (looks_like_number($TSD_size)) {
                            $TSD_type = "";
                            $entry_accepted=1;
                        }
                        else {
                            warn "WARNING: Your TSD coordinate entry is not valid\n";
                        }
                    }
                    else {
                        warn "WARNING: Your coordinate entries don't seem to be numbers\n";
                    }
                    unless ($TIR_size) {
                        $TIR_size = "N/A";
                    } 
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

            my $consensus_sequence; # this holds the consensus sequence and will be used to display the TIR sequence
            if ($menu1) { # only continue if the user has not elected to quit menu 1
                # The TIRs and TSDs location have now been selected, next create an alignment to display these to the user  
                my $temp_alignment_file = "/tmp/$element_name.fa"; # tried using perl temporary file system, but aliview will not open those
                open (OUTPUT, ">", $temp_alignment_file) or die "Cannot create temporary alignment file $temp_alignment_file\n";

                foreach my $seq_name (keys %alignment_sequences) {
                    if ($seq_name =~ /consensus/) { 
                        $consensus_sequence = $alignment_sequences{$seq_name};
                    }
                    else {# avoid the line with the consensus sequence
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
                `pkill java`; # kill a previous aliview window, this could be dangerous in the long run
                `aliview $temp_alignment_file`;
                if ($?) { die "Error executing: aliview $temp_alignment_file, error code $?\n"}

                my $menu2 = 1; # boolean, set to 1 until the user is done with menu 2
                my $element_rejected = 0; # boolean, set to 0 unless option "this is not an element selected", used to know which README to edit

                my $TIR1_sequence = substr($consensus_sequence, ($TIR_b1-1), $TIR_size);
                my $TIR2_sequence = substr($consensus_sequence, ($TIR_b2-$TIR_size), $TIR_size);
                while ($menu2) { #keep displaying until the user ready to leave
                    print color('green');
                    print "\nMENU 2 Select what to do with this element:\n";
                    print color('reset');
                    print "0) Go back to the previous menu\n";
                    # figure out the sequences of the current TIRs
                    
                    if ($TSD_type eq "TA") {
                        print "1) Update the README to say this is an element with TSDs of type TA and TIRs $TIR1_sequence and $TIR2_sequence\n";
                    }
                    else {
                        print "1) Update the README to say this as an element with TSDs of size $TSD_size and TIRs $TIR1_sequence and $TIR2_sequence\n";
                    }
                    print "2) Update the README to say this is not an element\n";
                    print "3) Change the TSD size\n";
                    print "4) Change the TIR sequence(s)\n";
                    print "5) Make a note in the README file\n";
                    print "6) Done reviewing this element\n";

                    do { # read the user input until it's a number within range
                        print "Line selection: ";
                        $pkey = <STDIN>;
                    } until ((looks_like_number($pkey)) and ($pkey <= 6));

                    # process the user's choice
                    if ($pkey == 0) {
                        $menu2 = 0;
                    } 
                    elsif ($pkey == 1) {
                        my $datestring = localtime(); 
                        if ($TSD_type eq "TA") {
                            print README "$datestring, Manual Review 1 result: This is an element, TSD TA, TIRs $TIR1_sequence and $TIR2_sequence\n";
                        }
                        else {
                            print README "$datestring, Manual Review 1 result: This is an element, TSD $TSD_size, TIRs $TIR1_sequence and $TIR2_sequence\n";
                        }
                    }
                    elsif ($pkey == 2) {
                        my $datestring = localtime(); 
                        print README "$datestring, Manual Review 1 result: This is not an element\n";
                        print REJECT "$datestring\t$element_name\tSTEP 3\tManual review of TSD and TIRs determined this is not an element\n";
                        `mv $ELEMENT_FOLDER/$element_name $reject_folder_path`;
                        if ($?) { die "ERROR: Could not move folder $ELEMENT_FOLDER/$element_name to $ANALYSIS_FOLDER: error code $?\n"}
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
                            else {
                                $TSD_type = ""; 
                            }
                        } until (looks_like_number($pkey2));
                        $TSD_size = $pkey2;
                    }
                    elsif ($pkey == 4) {
                        my $pkey2;
                        print "Enter TIR1 sequence [$TIR1_sequence]: ";
                        $pkey2 = <STDIN>;
                        chomp $pkey2;
                        if($pkey2) { # a new sequence has been entered
                            $TIR1_sequence = $pkey2;
                        }
                        print "Enter TIR2 sequence [$TIR2_sequence]: ";
                        $pkey2 = <STDIN>;
                        chomp $pkey2;
                        if($pkey2) { # a new sequence has been entered
                            $TIR2_sequence = $pkey2;
                        }
                    }
                    elsif ($pkey == 5) {
                        my $datestring = localtime(); 
                        if ($element_rejected) {
                            print "Edit the file $ANALYSIS_FOLDER/$element_name/README.txt starting with\n";
                        }
                        else {
                            print "Edit the file $ELEMENT_FOLDER/$element_name/README.txt starting with\n";
                        }
                        print "$datestring, Manual Review 1 user note: \n";
                        print "Press enter when done: ";
                        <STDIN>;
                    }
                    elsif ($pkey == 6) {
                        $menu2 = 0;
                        $menu1 = 0;
                    }
                }
            }            
        }
        close README;
    }
}

### PIPELINE STEP 4 
### Identify most likely transposase ORF based on identified TIR sequences

my $step_number = 4;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## Constant for this step
    my $MAX_ELEMENT_SIZE = 5000; # maximum element size
    
     ## update the analysis file with what is going on
    print ANALYSIS "Running STEP 4\n";
    unless ($INPUT_GENOME and $INPUT_PROTEIN_SEQUENCES){
        die "ERROR: for this step you need to provide two pieces of information:
             1) a fasta formated genome file, using the -g parameter
             2) the input protein file using the -p parameter\n";
    }
    print ANALYSIS "\tGenome: $INPUT_GENOME\n";
    print ANALYSIS "\tInput protein file: $INPUT_PROTEIN_SEQUENCES\n";
    print ANALYSIS "\tMAX_ELEMENT_SIZE = $MAX_ELEMENT_SIZE\n";

    ## Read the README files and identify TIRs sequences and TSD
    my %file_tirs;  # holds the element name as key and information on the element ends as array of values. 
                    # specifically [0] = TSD size or type, [1] = TIR1 sequence, [2] = TIR2 sequence
    opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
    while (readdir $dh) {
        unless ($_ =~ /^\./) {
            if (-e "$ELEMENT_FOLDER/$_/README.txt") {
                open (INPUT, "$ELEMENT_FOLDER/$_/README.txt") or die "ERROR: Cannot open file $!";
                while (my $line = <INPUT>) {
                    if ($line =~ /This\sis\san\selement,\sTSD\s(\S+),\sTIRs\s(\S+)\sand\s(\S+)/) {
                        $file_tirs{$_}[0] = $1; # TSD size and type
                        $file_tirs{$_}[1] = $2; # first TIR sequence
                        $file_tirs{$_}[2] = $3; # second TIR sequence
                    }
                }
                close INPUT; 
            }
            else {
                warn "WARNING: No README.txt file was found in folder $ELEMENT_FOLDER/$_";
            }
            unless (exists $file_tirs{$_}) {
                warn "WARNING: No TIR sequences were reported in the README file for element $_\n";
            }
        }
    }
    unless (keys %file_tirs) {
        warn "WARNING: No files were found that needed ORFs identified\n";
    }

    ## load all the protein sequences and change the name to match the folder names
    my %proteins = fastatohash($INPUT_PROTEIN_SEQUENCES);
    foreach my $name (keys %proteins) {
        if ($name =~ /^(\S+)/) {
            $proteins{$1}=$proteins{$name};
            delete $proteins{$name};
        }
    }

    ## Identify the position in the genome of all the TIR sequences for each element
    my %genome = fastatohash($INPUT_GENOME);
    foreach my $element_name (keys %file_tirs) { # go through the elements
       
        # identify the nucleotide sequences of sequences in between TIR locations
        my %element_sequences; # holds the genomic sequence of all the elements with intact tirs, it's a hash to avoid duplications
        foreach my $chr (keys %genome) { # go through each genome subsections (calling it "chr" here)
            # identify all the element sequences in the forward orientation. Converting everything to lower case to avoid confusion with cases.
            # Also provinding the name of the chrososome and orientation so that the position of all the elements can recorded
            my $tir1_seq = lc($file_tirs{$element_name}[1]);
            my $tir2_seq = lc($file_tirs{$element_name}[2]);
            %element_sequences = identify_element_sequence(lc($genome{$chr}), $tir1_seq, $tir2_seq, $MAX_ELEMENT_SIZE, $chr, "+");
            unless ($tir1_seq eq (rc($tir2_seq))) { # only look on the other strand if the TIRs are not symetrical, symetrical TIR will already have been found
                my %rc_element_sequences = identify_element_sequence(lc($genome{$chr}), lc($file_tirs{$element_name}[1]), lc($file_tirs{$element_name}[2]), $MAX_ELEMENT_SIZE, $chr, "-");
                foreach my $rc_element (keys %rc_element_sequences) {
                    $element_sequences{$rc_element} = $rc_element_sequences{$rc_element};
                }
            }

            foreach my $k (keys %element_sequences) {
                print "$k\n";
            }
#             for my $rce (keys %rc_element_sequences) {
# #                print "$rce\n";
#                 if ($rce =~ /(\S+):(\d+)-(\d+)/) {
#                     my $test_location = "$1:$3-$2";
#                     print "$test_location, $file_tirs{$element_name}[1], $file_tirs{$element_name}[2]\n";
#                     if (exists $element_sequences{$test_location}) {
#                         print "found one\n";
#                     }
                    
#                 }
  #          }
 #           my @reverse_sequence = identify_element_sequence(rc(lc($genome{$chr})), lc($file_tirs{$element_name}[1]), lc($file_tirs{$element_name}[2]), $MAX_ELEMENT_SIZE);
            # foreach my $sequence (@forward_sequence) {
            #     $element_sequences{$sequence} = 0;
            # }
            # foreach my $sequence (@reverse_sequence) {
            #     $element_sequences{$sequence} = 0;
            # } 
        }   

        # # make blast database file that contains the nucleotide sequences of the elements that were found
        # my $database_input_file = File::Temp->new(UNLINK => 1); # file that contains the nucleotide sequence of all the identified elements
        # open (OUTPUT, ">", $database_input_file) or die "$!";
        # foreach my $seq (keys %element_sequences) {
        #     print OUTPUT ">seq", rand(100000), "\n$seq\n"; # fasta lines with random name
        # }   
        # close OUTPUT;
        # my $tblastn_database_name = File::Temp->new(UNLINK => 1); # file name for the tblastn database name
        # `makeblastdb -in $database_input_file -dbtype nucl -out $tblastn_database_name`;
        # if ($?) { die "ERROR executing makeblastdb: error code $?\n"}

        # # make file with original protein used to find the the current element
        # my $protein_file = File::Temp->new(UNLINK => 1);
        # open (OUTPUT, ">", $protein_file) or die "$!";
        # print OUTPUT ">protein\n$proteins{$element_name}\n";
        # close OUTPUT;

        # # run tblastn
        # my $tblastn_output = File::Temp->new(UNLINK => 1); 
        # `tblastn -query $protein_file -db $tblastn_database_name -out $tblastn_output`;
        # if ($?) { die "ERROR executing makeblastdb: error code $?\n"}
        # `cp $tblastn_output ~/Desktop`;
        # print "done with $element_name\n";
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
    my ($chr_seq, $tir1, $tir2, $maximum_size, $chromosome_name, $orientation) = @_; # takes as input the chromosome sequence, the TIR sequences,
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
    # identify any pairs that are in the correct orientation and not too far from on another
    for (my $i=0; $i<scalar @tir1; $i++) {
        for (my $j=0; $j<scalar @tir2; $j++) {
            my $size = $tir2[$j] - $tir1[$i];
            if (($tir1[$i] < $tir2[$j]) and ($size < $maximum_size)){
                my $b1; # boundary of the element on the chromosome
                my $b2; # boundary of the element on the chromosome
                if ($orientation eq "+") {
                    $b1 = $tir1[$i]+1; 
                    $b2 = $tir2[$j];
                }
                else {
                    $b2 = $tir1[$i]+1;
                    $b1 = $tir2[$j];
                }
                my $location_name = "$chromosome_name:$b1-$b2"; 
                $nucleotide_sequences{$location_name} = substr($chr_seq, $tir1[$i], $size);
            }
        }
    }
    return (%nucleotide_sequences);
}