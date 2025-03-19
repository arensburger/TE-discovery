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

1;