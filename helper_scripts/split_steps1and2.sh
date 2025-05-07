#!/bin/sh
# The purpose of this script is to run the pipeline through steps 1 and 2, but using as 
# input a large tblastn file (such as one created by running the following command on a 
# cluster:
#
# tblastn -query <proteins fasta> -db <genome fasta> -outfmt "6 qseqid sseqid sstart send pident length qlen" -out <tblastn.o> 
# (note, it's suggested that your first randomize the <protein fasta> file, using the "randomfasta.pl" script, to 
# increase the chances that the workload will be evenly distributed among sub files)

# RUNNING THIS SCRIPT
# This script creates a number of files, to avoid overwriting any of them it's best
# to create a new temporary directory in the same folder as this script.

# Update the variables below to suit the current analysis
TEMP_FOLDER="./temp" # path to the new temporary directory
TBLASTN_LOC="./tblastn.o" # path to the tblastn.o output file
NUM_SECTIONS=6 # number of sections to split the input file into
SCRIPT_LOC="/home/peter/TE-discovery/mainscript.pl" # location of mainscript.pl
GENOME_LOC="/home/peter/db/Ptep_3.0/GCF_000365465.3_Ptep_3.0_genomic.fna" # location of genome, formated with 'makeblastdb', and file .length has been generated
ANALYSIS_FOLDER="./Ptep3_Analysis_files" # path where to store the combined results
ELEMENT_FOLDER="./Ptep3_Element_files" # path where to store the combined results

# # STEP 1:	Create the temporary directory, randomize the tblasn.o output, nd split it in multiple files
#mkdir $TEMP_FOLDER # make a new folder
#split -n l/$NUM_SECTIONS $TBLASTN_LOC $TEMP_FOLDER/ # split the file

# # STEP 2:	Run the pipeline steps 1 and 2
# Uncoment the appropriate number of lines depending on how many sections you have
#perl $SCRIPT_LOC -t $TEMP_FOLDER/aa -g $GENOME_LOC -n $TEMP_FOLDER/analysis1 -e $TEMP_FOLDER/element1 -a 1 -b 2 &
#perl $SCRIPT_LOC -t $TEMP_FOLDER/ab -g $GENOME_LOC -n $TEMP_FOLDER/analysis2 -e $TEMP_FOLDER/element2 -a 1 -b 2 &
#perl $SCRIPT_LOC -t $TEMP_FOLDER/ac -g $GENOME_LOC -n $TEMP_FOLDER/analysis3 -e $TEMP_FOLDER/element3 -a 1 -b 2 &
#perl $SCRIPT_LOC -t $TEMP_FOLDER/ad -g $GENOME_LOC -n $TEMP_FOLDER/analysis4 -e $TEMP_FOLDER/element4 -a 1 -b 2 &
#perl $SCRIPT_LOC -t $TEMP_FOLDER/ae -g $GENOME_LOC -n $TEMP_FOLDER/analysis5 -e $TEMP_FOLDER/element5 -a 1 -b 2 &
#perl $SCRIPT_LOC -t $TEMP_FOLDER/af -g $GENOME_LOC -n $TEMP_FOLDER/analysis6 -e $TEMP_FOLDER/element6 -a 1 -b 2 &

# STEP 3:	Combine the results into element and analysis folders (run this after STEP 2 is done)
 mkdir $ANALYSIS_FOLDER
 mkdir $ELEMENT_FOLDER
 mkdir "$ANALYSIS_FOLDER/Rejected_elements"

 i=1
 while [ $i -le $NUM_SECTIONS ]
 do 
	echo "Working on $i"
 	cp -r $TEMP_FOLDER/analysis$i/Rejected_elements/* "$ANALYSIS_FOLDER/Rejected_elements"	
 	cat $TEMP_FOLDER/analysis$i/Analysis_parameters.txt >> $ANALYSIS_FOLDER/Analysis_parameters.txt
 	cat $TEMP_FOLDER/analysis$i/Rejected_sequences.txt >> $ANALYSIS_FOLDER/Rejected_sequences.txt
 	cp -r $TEMP_FOLDER/element$i/* $ELEMENT_FOLDER	
 	i=$(($i+1))
done

# echo "Done. Temporary directory $TEMP_FOLDER can be deleted. Results have been copied to $ANALYSIS_FOLDER and $ELEMENT_FOLDER"
