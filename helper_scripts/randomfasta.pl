#!/usr/bin/perl

# June 2015 samples random sequences from a fasta formated file
# May 2025 updating

use strict;
use File::Temp ();
use Getopt::Long;
use List::Util 'shuffle';

##> Define options
my %config;
GetOptions (\%config,
            'in=s',
            'n=s',
	    	'o=s');

##> Check if no mandatory parameter is missing and set defaults
if (!exists $config{in}) {printUsage();}
if (!exists $config{o}) {printUsage();}

##> load the data and randomize 
my %sequences = genometohash($config{in}); # load the input file
my @fasta_titles = keys %sequences; # get the fasta titles (assumes they are unique)
my @randtitle = shuffle @fasta_titles; # randomize array

##> Check on the optional input about the number of sequences to output
my $numseq = scalar @fasta_titles;
if ((exists $config{n}) and ($config{n} > $numseq )) {
	die "ERROR: The number of seqeuences to output ($config{n}) is bigger than the number of sequences in the file ($numseq)\n";
}
elsif (!exists $config{n}) {
	$config{n} = $numseq;
}

##> Print output
open (OUTPUT1, ">$config{o}") or die "cannot open output file $config{o}\n";
for (my $i=0; $i<$config{n}; $i++) {
	print OUTPUT1 ">", $randtitle[$i], "\n";
	print OUTPUT1 $sequences{$randtitle[$i]}, "\n";
}
close OUTPUT1;

###########################################################################
################################ Functions ################################
###########################################################################


###########################################################################
sub printUsage{

print STDOUT "USAGE : perl randomfasta.pl -in \"fasta file\" -n \"number of sequences to sample\" -o \"output file\"
Options : 
    -in 	fasta file where each sequence is only one line long (Mandatory)
    -n 		number of sequences to sample, if no number provided return all the sequences randomised (Optional)
    -o    	Name of ouput file (Mandatory)\n";
    exit;
}

###########################################################################
#load a genome into a hash
sub genometohash {
	use strict;
	(my $filename) = @_;
	my %genome; #hash with the genome
	my $seq="";
	my $title;
	open (INPUT100, $filename) or die "cannot open input file $filename in sub genometohash\n";
	my $line = <INPUT100>;
	my $title;
	my $seq = "";

	# dealing with the first line
	if ($line =~ />(.+)/) {
		chomp $1;
		$title = $1;
	}
	else {
		die "ERROR, fasta file $filename does not start with a title line\n";
	}
	## dealing with the remaining lines 
	while (my $line = <INPUT100>) {
		if ($line =~ />(.+)/)  {
			chomp $1;
			my $new_title = $1;
			if (exists $genome{$new_title}) {
				print STDERR "error in sub genometohash, two contigs have the name $new_title, ignoring one copy\n";
			}
			else {
				$genome{$title} = $seq;
				$title = $new_title;
			}
			$seq = "";
		}
		else {
			$line =~ s/\s//g;
			$seq .= $line;
		}
	}
	$genome{$title} = $seq;
	close INPUT100;
	return (%genome);
}