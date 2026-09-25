# September 2026. Convert cluster report files to RepeatMasker format

use strict;
use Getopt::Long;
use File::Find;

my $CLUSTER_FOLDER; # folder with the reports
### GET INPUTS Read and check that the inputs have been provided
GetOptions(
	'c:s'   => \$CLUSTER_FOLDER
);

## CHECK INPUTS Validate required input files
unless ($CLUSTER_FOLDER) {
    die "usage: perl cluster2rm.pl <-c path to folder with cluster report files REQUIRED>\n";
}
$CLUSTER_FOLDER = fixdirname($CLUSTER_FOLDER);

## Identify all the report files in the cluster
my @reports; # holds the full path name to the report files
find(sub {
    return unless -f $_;
    if ($_ =~ /^cluster\d+_report\.md$/) {
        push @reports, $File::Find::name;
    }
}, $CLUSTER_FOLDER);

## Check that at least one report was found
unless (@reports) { die "ERROR: No cluster reports were found in the path provided\n"}

## Read the reports and outuput RepeatMasker formated output
foreach my $filename (@reports) {
    my $cluster_name;
    my $nucleotide_sequence; 
    my $taxonomy;
    open (INPUT, $filename) or die "ERROR: Cannot open file $filename $!\n";
    while (my $line = <INPUT>) {
        if ($line =~ /# Report (cluster\d+)/) {
            $cluster_name = $1;
        }
        if ($line =~ /^## Nucleotide sequence:/) {
            $nucleotide_sequence = <INPUT>;
            chomp $nucleotide_sequence;
        }
        if ($line =~ /^## Element taxonomy:/) {
            $taxonomy = <INPUT>;
            chomp $taxonomy;
        }
    }

    # check that all the fields were found, report an error if not
    if ($cluster_name and $nucleotide_sequence and $taxonomy) {
        print ">$cluster_name#DNA/$taxonomy\n$nucleotide_sequence\n";
    }
    else {
        warn ("The file $filename seems to be missing information, no RepeatMasker output was created\n")
    }
    
}

# Takes a folder name and standardizes the formating
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