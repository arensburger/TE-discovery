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
    # if ($min_tir_size <= ($pos -1)) {
    #     my $tir1_sequence = substr($s1, 0, ($pos - 1));
    #     my $tir2_sequence = substr($s2, -($pos - 1), ($pos - 1));
    #     return ($tir1_sequence, $tir2_sequence);
    if ($min_tir_size <= $lastgoodbase) {
        my $tir1_sequence = substr($s1, 0, ($lastgoodbase+1));
        my $tir2_sequence = substr($s2, -($lastgoodbase+1), ($lastgoodbase+1));
        return ($tir1_sequence, $tir2_sequence);
    }
    else {
        return ("","");
    }

	# ### get the ends into string Variables
	# my $s1 = substr ($sequence, 0, $SCAN_SIZE);
	# my $s2 = substr ($sequence, -$SCAN_SIZE, $SCAN_SIZE);
	# # reverse complement seq2
	# my $s2rc = rc($s2);
	# ### count the matches
	# my $matches;
	# for (my $i=0; $i < length $s1; $i++){
	# 	my $char_s1 = substr ($s1, $i, 1);
	# 	my $char_s2 = substr ($s2rc, $i, 1);

	# 	my $countN_char1 = () = $char_s1 =~ /N|n|-/; # count forbiden characters
	# 	my $countN_char2 = () = $char_s2 =~ /N|n|-/;
	# 	if (($i==0) and ($countN_char1 or $countN_char2)) { # check if the TIR starts with forbidden character
	# 		return (0, "", "");
	# 	}

	# 	unless ($countN_char1 or $countN_char2) {
	# 		if ($char_s1 eq $char_s2) {
	# 			$matches++;
	# 		}
	# 	}
	# }
	# ### evaluate the results
	# if ($matches >= $MIN_MATCH) {
    #     my $c1 = substr($s1, 0, 3);
    #     my $c2 = substr($s2, -3, 3);
	# 	return (1, $c1, $c2); # tir found, report that it was found and bases
	# }
	# else {
	# 	return (0, "", ""); # tir not found
	# }
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