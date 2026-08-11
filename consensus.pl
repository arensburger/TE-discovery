#!/usr/bin/ perl
use strict;
use Getopt::Long;

### CONSTANTS
my $CONSENSUS_LEVEL = 0.6; # minimum proportion of non-gap positions that must agree to call a position
my $MIN_FOR_POSITION = 2; # minumum number of nucleotides or amino acids to call a position (not ignore it)
my $TIRTSD; # optional, report the TSDs and TIRs based on the sequence names
my $SHOW_HELP; # call for help 

### Set up input parameters
GetOptions(
	'c:s'   => \$CONSENSUS_LEVEL,
    'p:s'   => \$MIN_FOR_POSITION,
    't'   => \$TIRTSD,
    'h'     => \$SHOW_HELP,
);
if ($SHOW_HELP) {
    print STDERR "This script takes aligned fasta formatted sequences in the clipboard as input and produces a consensus sequence. Optionally it can ouptut TSDs and TIR sequences as well. If the TSD and TIR sequences are requested, the consensus ouput sequence does not include TSD sequences.\n";
    print STDERR "usage: perl consensus.pl <-t display TIR and TSD information based on sequence name OPTIONAL> <-c set consensus level, default $CONSENSUS_LEVEL OPTIONAL> <-p set mininum number of sequences for consensus, default $MIN_FOR_POSITION OPTIONAL>\n";
    exit;
}


### Read the clipboard for fasta formatted sequences and verify that all is ok
my $text = `xclip -selection clipboard -o 2>/dev/null`;
if (!defined $text || $text eq '') {
    die "Could not read clipboard. Make sure 'xclip' is installed.\n";
}
die "Clipboard is empty or could not be read.\n" unless defined $text && length $text;

### Parse the fasta sequence
my %sequences;
my $header;
for my $line (split /\r?\n/, $text) {
    if ($line =~ /^(>\S+)/) {
        $header = $1;
        $sequences{$header} = '';
    }
    elsif (defined $header && $line =~ /\S/) {
        $sequences{$header} .= $line;
    }
}

### Check the content of %sequences
if (!%sequences) { # check if FASTA sequences were read
    die "ERROR: No FASTA formated sequences were found in the clipboard.\n";
}
elsif (scalar keys %sequences == 1) { # check if there's more than one sequence
    my ($header) = keys %sequences;
    die "ERROR: The clipboard only appears to have a single FASTA formatted entry, here's what was read\n$header\n$sequences{$header}\n";
}
## check if the sequences all have the same length
my @lengths = map { length $_ } values %sequences;
my %unique_lengths = map { $_ => 1 } @lengths;
unless (scalar keys %unique_lengths == 1) {
    die  "ERROR: The input FASTA sequences don't all have the same length, need to have aligned sequences as input\n";
} 

### figure out if these are nucleotides or proteins 
my $single_sequence; 
foreach my $s (keys %sequences) {
    $single_sequence .= $sequences{$s};
}
my $type = guess_seq_type($single_sequence);

### Create consensus sequence
my $consensus; 
for (my $i=0; $i < $lengths[0]; $i++ ) { # go through each position one at a time
    
    ## Store the characters at this postion into an arry
    my @chars; # array that holds all the characters at this position
    for my $header (keys %sequences) { # go through each sequence at this position
        my $c = substr $sequences{$header}, $i, 1; # current character
        unless ($c eq "-") { # don't add gap characters
            push @chars, $c;
        }
    } 

    ## Identify the most abundant character
    my %count;
    $count{$_}++ for @chars;
    my ($top_char) = sort { $count{$b} <=> $count{$a} } keys %count;
    my $total_characters = scalar @chars;
    my $proportion = $count{$top_char} / $total_characters;

    ## Decide if this character should be printed or not
    if(($total_characters >= $MIN_FOR_POSITION) and ($proportion >= $CONSENSUS_LEVEL)) { # condiditions to report a non N consensus
        $consensus .= "$top_char";
    }
    elsif (($total_characters >= $MIN_FOR_POSITION) and ($proportion < $CONSENSUS_LEVEL)) { # conditions to report N or x
        if ($type eq 'nucleotide') {
            $consensus .= "N";
        }
        else {
            $consensus .= "X";
        }
    }
}

### Print output
my $TSDseq;
my $TSDlength;
my $TIR1seq; 
my $TIR2seq;
if ($TIRTSD and $consensus) { # user has requested the TIR and TSD information be displayed, and a consensus must have been created
    my ($header) = keys %sequences; # random header with the TSD and TIR information   
    if ($header =~ /\d+_(\d+)_(\d+)$/) {
        $TSDseq = uc(substr $consensus, 0, $1);
        $TIR1seq = uc(substr $consensus, $1, $2);
        $TIR2seq = uc(substr $consensus, -$1-$2, $2);     
        if ($TSDseq eq "TA") { 
            print "TSD: TA\n";
            $TSDlength = 2;
        } 
        else {
            $TSDlength = length $TSDseq;
            print "TSD: $TSDseq.bp\n";
        }
        print "TIR1: $TIR1seq\n";
        print "TIR2: $TIR2seq\n";
    }
    else {
        die "WARNING: Could not extract TSD and TIR information from header: $header\n";
    }
}

if ($consensus) { # print the consensus sequence if the TSD and TIR was requested then don't print the TSDs
    my $ali_length = $lengths[0] - (2*$TSDlength);
    my $cons_length = (length $consensus) - (2*$TSDlength);
    my $num_input_seq = keys %sequences;
    print "Alignment length: $ali_length, Consensus length $cons_length, Number of input sequences: $num_input_seq\n";
    print ">consensus\n";
    if ($TIRTSD) { 
        my $no_tsd_cons = substr $consensus, $TSDlength, (length $consensus) - (2*$TSDlength);
        print "$no_tsd_cons\n";
    }
    else {
        print "$consensus\n";
    }
}
else {
    print "No consensus found for these input sequences:\n";
    for my $header (keys %sequences) {
        print "$header\n$sequences{$header}\n";
    }
}

# figure out type
sub guess_seq_type {
    my ($seq) = @_;
    $seq = uc $seq;
    $seq =~ s/[^A-Z]//g;  # strip non-letters (gaps, whitespace, etc.)

    my $len = length $seq;
    return 'unknown' unless $len;

    my $nt_chars = ($seq =~ tr/ACGTNacgtn//);
    my $nt_fraction = $nt_chars / $len;

    # Nucleotide sequences are almost entirely A/C/G/T/N
    return $nt_fraction >= 0.9 ? 'nucleotide' : 'protein';
}