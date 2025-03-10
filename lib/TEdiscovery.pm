### These are scripts required for TE-discovery 

### Takes a fasta file name as input and return a hash with content
sub fastatohash {
    (my $filename) = @_;
    my %seqhash; #final hash with the seqhash
    my $current_header; # fasta header of the currently read sequence

    open (INPUT, $filename) or die "ERROR: cannot open input file $filename in fastatohash subroutine\n";
    my $line = <INPUT>;

    ## record the first header
    if ($line =~ /^>(.+)/) {
        $current_header = $1;
    }
    else {
        die "ERROR: file $filename does not appear to be a FASTA formatted file, this in the fastatohash subroutine\n";
    }

    while ($line = <INPUT>) {
        if ($line =~ /^>(.+)/) {
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

1;