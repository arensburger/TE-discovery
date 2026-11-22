# Octber 2026. Converts cluster report files to just the TIR sequences

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

## Read the reports a
foreach my $filename (@reports) {
    my $TIR1;
    my $TIR2;
    my $cluster_name; 
    open (INPUT, $filename) or die "ERROR: Cannot open file $filename $!\n";
    while (my $line = <INPUT>) {
        if ($line =~ /#\sReport\s(\S+)/) {
            $cluster_name = $1;
        }
        if ($line =~ /-\s5'\sTIR:\s(\S+)/) {
            $TIR1 = $1;
        }
        if ($line =~ /-\s3'\sTIR:\s(\S+)/) {
            $TIR2 = $1;
        }
    }

    # check that all the fields were found, report an error if not
    if ($cluster_name and $TIR1 and $TIR2) {
        print ">$cluster_name-TIR1\n$TIR1\n";
        print ">$cluster_name-TIR2\n$TIR2\n";
    }
    else {
        warn ("The file $filename seems to be missing information\n")
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