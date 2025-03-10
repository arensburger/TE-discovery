# Feb 2025 - Top-level script for TE discovery pipeline
# Inputs: List of TE sequences (FASTA) and a genome (FASTA)

use strict;
use Getopt::Long;
use File::Temp qw(tempfile);
use FindBin::libs;  # this sets up that the directory lib will have all modules necessary to run the program 
use TEdiscovery;    # these are the subcripts necessary to run the pipeline
use List::UtilsBy qw(max_by);

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
my $GAP_THRESHOLD=0.75; # REMOVING GAPS FROM ALIGNMENT, if an alignment position has this proportion or more of gaps, then remove it from the multiple sequence alignment
my $CONSLEVEL=0.60; # MAKING CONSENSUS OF SEQUENCES sequence consensus level for consensus
my $WINDOW_SIZE = 15; # MAKING CONSENSUS OF SEQUENCES size of the window looking for stretches of N's
my $SCAN_SIZE = 100; # MAKING CONSENSUS OF SEQUENCES number of bp to scan each side of alignment to determine if there are N's on the edges at all
my $SCAN_PROP = 0.5; # MAKING CONSENSUS OF SEQUENCES mininum proportion of side scan that has to be N's to be a real N full edge
my $MAX_WIN_N = 2; # MAKING CONSENSUS OF SEQUENCES maximum number of N's in the first window where the transition from N to non-N is
my $EDGE_TEST_PROPORTION = 0.05; # TESTING CONSENSUS SEQUENCES how far from the edge of the consensus do the non-gap positions have to start

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
my $ANALYSIS_FILES_OUTPUT_DIR="./$ANALYSIS_NAME-analysis-files"; # directory to store output files of current analysis (no slash at the end), these can be destroyed when analysis is finished
if (-d $ANALYSIS_FILES_OUTPUT_DIR) {
    warn "PRELIMINARY STEPS: Directory $ANALYSIS_FILES_OUTPUT_DIR already exists, using it to store files generated during this analysis (existing files with the same name will be overwriten)\n";
}
else {
    warn "PRELIMINARY STEPS: Creating directory $ANALYSIS_FILES_OUTPUT_DIR for storing files generated during the analysis";
    mkdir( $ANALYSIS_FILES_OUTPUT_DIR ) or die "Couldn't create $ANALYSIS_FILES_OUTPUT_DIR directory, $!";
}

## CHECK INPUTS Create output directory for individual elements if necessary
my $ELEMENT_FOLDER="./$ANALYSIS_NAME-elements"; # directory where the analysis of individual element will be stored
if (-d $ELEMENT_FOLDER) {
    warn "PRELIMINARY STEPS: Directory $ELEMENT_FOLDER already exists, using it to put folders for individual elements (any current subfolders will not be overwritten)\n";
}
else {
    warn "PRELIMINARY STEPS: Creating directory $ELEMENT_FOLDER for that will have subdirectories for individual elements";
    mkdir( $ELEMENT_FOLDER ) or die "Couldn't create $ELEMENT_FOLDER directory, $!";
}

## VARIABLES used by more than one step in the pipeline
my $blast_output_file_name = "$ANALYSIS_FILES_OUTPUT_DIR/tblastn.o";

### PIPELINE STEP 1 identify proteins that match the genome with parameters specified above under "CONSTANTS"
###     The output is a list of proteins for further analysis recorded in the file $output_file_name
my $step_number = 1;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## VARIABLES, variable for this step
    my %protein_ids; # holds the id the input proteins that passed the filtering tests as key and the number of copies that passed the test as values

    ## Excute the tblastn search
   `tblastn -query $INPUT_PROTEIN_SEQUENCES -db $INPUT_GENOME -outfmt "6 qseqid sseqid sstart send pident length qlen" -out $blast_output_file_name -num_threads $NUM_THREADS`;
   if ($?) { die "ERROR executing tblastn, stoping analysis (hint: was the genome formated with makeblastdb?): error code $?\n"}

    ## Inspired by the Goubert et al. protocol, filter elements that 1) have >= 80% identity to genome, 2) have 50% length of the query, 3) are found at multiple locations
    my %candidate_protein; # hash with protein name as key and string with chromosome and middle location of element on that chromsome

    open (INPUT, "$blast_output_file_name") or die "ERROR: Cannot open file $blast_output_file_name\n";
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

    # Filter out elements that have too few copy numbers
    # (making this a separate step so the code is more modular, rather then incoroporating it into the next step)
    foreach my $prot_name (keys %protein_ids) {
        unless ($protein_ids{$prot_name} >= $COPY_NUMBER) {
            delete $protein_ids{$prot_name};
        }
    }

    # Create individual directories for each element
    my $i=0; # counts the number of output lines, to check if it's zero
    foreach my $prot_name (keys %protein_ids) {
        my $folder_name = "$ELEMENT_FOLDER/$prot_name";
        if (-d $folder_name) {
            warn "STEP $step_number: Element folder $folder_name already exists, using it without modification to existing files\n";
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
my $step_number = 2;
if (($step_number >= $START_STEP) and ( $step_number <= $END_STEP)) { # check if this step should be performed or not  
    print "Working on STEP $step_number ...\n";

    ## VARIABLES, variable for this step
    my @elements; # name of all the elements that will be anlaysed in this step

    ## STEP 2.1
    ## figure out the elements that we're working with 
    ## (for the sake of being modular, redoing this instead of just taking the data from the previous step)
    opendir(my $dh, $ELEMENT_FOLDER) or die "ERROR: Cannot open element folder $ELEMENT_FOLDER, $!";
    while (readdir $dh) {
        unless ($_ =~ /^\./) { # prevents reading invivisible files or the . and .. files
            push @elements, $_;
        }
    }
    unless (scalar @elements) { # check that at least one element is present to analyse
        die "ERROR: No elements to analyse were found in the folder $ELEMENT_FOLDER\n";
    }

    ## STEP 2.2
    ## identify TIRs
    my $i=1;
    foreach my $element_name (@elements) {

        #report progess to screen
        my $size = scalar @elements;
        print "\tProcessing $element_name, number $i of $size\n";
        $i++;

        # create or open the README file for this element
        open (README, ">$ELEMENT_FOLDER/$element_name/README.txt") or die "ERROR: Could not open or create README file $ELEMENT_FOLDER/$element_name/README.txt\n";

        # STEP 2.2.1
        # extend the blast hits, the output of this step is a fasta file ending in -extend.fa
        my @blastlines = (); # holds all the relevant blast lines for this element
        open (INPUT, $blast_output_file_name) or die "ERROR: cannot open file blast file output $blast_output_file_name\n";
        while (my $line = <INPUT>) {
            my @data = split ' ', $line;
            if ($data[0] eq $element_name) {
                push @blastlines, $line;
            }
        }
        close INPUT;

        # checking that the genome length file is present
        unless (-f "$INPUT_GENOME.length") {
            die "ERROR: Cannot file find file $INPUT_GENOME.length needed for bedtools, generate one using\n(assuming it's been formated for BLAST with makeblastdb)\n\nsamtools faidx \$genome\nawk \'{OFS=\"\\t\"; print \$1,\$2}\' < \$genome.fai > \$genome.length\n";
        }

        # create the bed file
        my $bedfile = File::Temp->new(UNLINK => 1, SUFFIX => '.bed' );
        open (OUTPUT, '>', $bedfile) or die "$!\n";
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

        # create the files with extended boundaries
        my $slopfile = File::Temp->new(UNLINK => 1, SUFFIX => '.slop' );
        `bedtools slop -s -i $bedfile -g "$INPUT_GENOME.length" -b $BLAST_EXTEND > $slopfile`;
        if ($?) { die "ERROR executing bedtools: error code $?\n"}

        `bedtools getfasta -fi $INPUT_GENOME -fo $ELEMENT_FOLDER/$element_name/$element_name-extend.fa -bed $slopfile -s`;
        if ($?) { die "ERROR executing bedtools: error code $?\n"}


        # # merge the overlaping elements into temporary file
        # `perl $SCRIPTS_FOLDER/process_folders-merge_overlaps2.pl -f $element-extend.fa -p $ELEMENT_FOLDER > temp-merged.fa`;
        # if ($?) { die "Error executing script $SCRIPTS_FOLDER/process_folders-merge_overlaps2.pl , error code $?\nNote: Check that the element folder names match the names in the blast file\n"}

        # STEP 2.2.2
        # align the sequences
        my $temp_aligned_sequences = File::Temp->new(UNLINK => 1, SUFFIX => '.maf' );
        `mafft --quiet --thread -$NUM_THREADS $ELEMENT_FOLDER/$element_name/$element_name-extend.fa > $temp_aligned_sequences`;
        if ($?) { die "Error executing mafft, error code $?\n"}

        # STEP 2.2.3
        # remove positions from the genome that have more gaps than the threshold at $GAP_THRESHOLD    
        my %cp; # cp = current position, sequence name as key and current position as value at [0] and orientation at [1]
                        # [2] is boolean, value 0 until a non-gap has been reached
        my %seq; # sequence name as key and string with mafft alignment as value
        my %seqrmg; # sequences with gaps removed
        my %pos; # holds the sequence name as key array with position of every nucleotide as value
        my $alilen; # length in bp of the multiple sequence alignment
        my $ali_trimmed_length; # length in bp of the multiple squence alignment after it's been trimmed

        # populate the %cp hash with information about the multiple sequence alignment
#$temp_aligned_sequences = "/home/peter/Desktop/test.maf";
        my %seq = fastatohash($temp_aligned_sequences); 
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
        close OUTPUT;

        # output the file that will allow restoring the deleted positions, called .alipos
        open (OUTPUT, '>', "$ELEMENT_FOLDER/$element_name/$element_name.alipos") or die "Error: cannot create file $element_name.alipos, $!\n";
        for my $name (keys %pos) {
            for my $j ( 0 .. $#{ $pos{$name}}) {
                my $location = $j+1;
                print OUTPUT "$name\t$location\t$pos{$name}[$j]\n";
            }
	    }
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

        # Test to see if this consensus sequences passes the "edge test". The trimmed consensus should not too close 
        # to the edge of the alignment (this would suggest that the whole region, not just an element inside it, is conserved)
        my $edge_test; # boolean, will be set to 1 if the consensus is not too close to the edge
        if (($ltrans > ($EDGE_TEST_PROPORTION * (length $trimmed_conseq))) and ($rtrans > ((1 - $EDGE_TEST_PROPORTION) * (length $trimmed_conseq)))) {
            $edge_test = 1;
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


open (OUTPUT, '>', "/home/peter/Desktop/ali.maf") or die "$!";
print OUTPUT ">cons\n";
print OUTPUT "$conseq\n";
print OUTPUT ">cons-$ltrans-$rtrans\n";
print OUTPUT "$trimmed_conseq\n";
foreach my $key (keys %seqrmg) {
    print OUTPUT ">$key\n";
    print OUTPUT "$seqrmg{$key}\n";
}
close OUTPUT;

        # STEP 2.2.6 
        # identify the TIR and TSD locations
        if ($trimmed_conseq and $edge_test) { # only continue if consensus was created and edge test was passed
            print "identify TIRs for $element_name\n"
        }
        else {
            print "edge test: $edge_test\n$trimmed_conseq\n";
        }

        # # make consensus and record boundaries of consensus sequence
        # `perl $SCRIPTS_FOLDER/make_consensus2.pl -in temp-aligned.maf2 > temp.cons`;
        # if ($?) { die "Error executing script $SCRIPTS_FOLDER/make_consensus.pl -in temp-aligned.maf2, error code $?\n"}

        # # test to see if consensus was created and if so load the name and sequence into memory
        # my $cons_b1; # left boundary position of consensus
        # my $cons_b2; # right boundary position of consensus
        # my $cons_sequence; # sequence of the consensus
        # if (-e "temp.cons") { # consensus exists
        #     open (INPUT, "temp.cons") or die "Cannot open file temp.cons in current directory\n";
        #     my $line = <INPUT>;
        #     if ($line =~ /^>consensus_(\d+)_(\d+)/) {
        #         $cons_b1 = $1;
        #         $cons_b2 = $2;
        #     }
        #     else {
        #         die "Error consensus file temp.cons has the title below that cannot be parsed\n$line";
        #     }
        #     my $line = <INPUT>;
        #     chomp $line;
        #     $cons_sequence = $line;
        # }
        # else { # consensus absent
        #     my $datestring = localtime();
        #     print OUTPUT "$datestring, ran script process_folders-blast-and-cons2.pl --> consensus sequence was NOT created\n";

        # }

        # # merge the consensus with the alignment
        # `cat temp.cons temp-aligned.maf2 > temp-aligned.maf3`;
        # if ($?) { die "Error executing script cat temp.cons temp-aligned.maf2 > temp-aligned.maf3, error code $?\n"}

        # # test if the edge test is passed, if it is then continue the analysis
        # if (($cons_b1 > ($LEFT_EDGE_PROP * (length $cons_sequence))) and ($cons_b2 < ($RIGHT_EDGE_PROP * (length $cons_sequence)))) {
        #     # update README with edge test results
        #     my $datestring = localtime();
        #     print OUTPUT "$datestring, ran script process_folders-blast-and-cons2.pl --> consensus edge test passed\n";

        #     # find tirs and tsds
        #     `perl $SCRIPTS_FOLDER/find_tir3.pl -a temp-aligned.maf3 -b $cons_b1 -c $cons_b2 -d 10 > temp.tir`;
        #     if ($?) { die "Error executing script perl $SCRIPTS_FOLDER/find_tir3.pl -a temp-aligned.maf3, error code $?\n"}

        #     `perl $SCRIPTS_FOLDER/process_folders-sort-tsds2.pl -f $element_name -t temp.tir -p $ELEMENT_FOLDER >> $EDGETESTPASS_FILE`;
        #     if ($?) { die "Error executing script perl $SCRIPTS_FOLDER/process_folders-sort-tsds2.pl -f $element_name -t temp.tir -p $ELEMENT_FOLDER, error code $?\n"}

        #     `mv temp.tir $ELEMENT_FOLDER/$element_name/$element_name.tir`;
        #     if ($?) { warn "Error code when moving file: mv $ELEMENT_FOLDER/$element_name/$element_name.tir $ELEMENT_FOLDER/$element_name, $?\n" }

        #     print OUTPUT "\tfile $element_name.tir: result of looking for tirs and tsds\n";
        # }
        # else { # if gets here then failed the consensus edge test
        #     push @no_pass, $element_name;

        #     # update README with edge test results
        #     my $datestring = localtime();
        #     print OUTPUT "$datestring, ran script process_folders-blast-and-cons2.pl  --> consensus edge test failed\n";
        # }

        # # print the rest of the README file
        # print OUTPUT "\tfile $element_name-extend.fa: result of blast file $blastfile extended by $BLAST_EXTEND bp\n";
        # print OUTPUT "\tfile $element_name.alipos: file with information to reconstitute aligment after removing 80% gaps\n";
        # print OUTPUT "\tfile $element_name.maf: blast alignment file with 80% gaps removed and consensus file added\n";
        close README;

        # # clean up
        # `mv temp-aligned.maf3 $ELEMENT_FOLDER/$element_name/$element_name.maf`;
        # if ($?) { warn "Error code when moving file:mv temp-aligned.maf3 $ELEMENT_FOLDER/$element_name/$element_name.maf, $?\n" }

        # `mv $element_name-extend.fa $ELEMENT_FOLDER/$element_name`;
        # if ($?) { warn "Error code when moving file: mv $element_name-extend.fa $ELEMENT_FOLDER/$element_name, $?\n" }

        # `rm temp-merged.fa`;
        # if ($?) { warn "Error code when removing file: rm temp-merged.fa, $?\n" }

        # `rm temp-aligned.maf`;
        # if ($?) { warn "Error code when removing file: rm temp-aligned.maf, $?\n" }

        # `rm temp.cons`;
        # if ($?) { warn "Error code when removing file: rm temp.cons, $?\n" }

        # `rm temp-aligned.maf2`;
        # if ($?) { warn "Error code when removing file: rm temp-aligned.maf2, $?\n" }

        # `mv temp-aligned.alipos $ELEMENT_FOLDER/$element_name/$element_name.alipos`;
        # if ($?) { warn "Error code when moving file: mv temp-aligned.alipos $ELEMENT_FOLDER/$element_name/$element_name.alipos, $?\n" }

    }    
}
