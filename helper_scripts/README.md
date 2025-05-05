- May 4, 2025 These scripts are to help run the pipeline
1) `randomfasta.pl` randomizes the order of the fasta file. This is useful to make sure that the files in the script `split_steps1and2.sh` will run more less at the same speed.
2) `split_steps1and2.sh` is designed to take a large tblastn output file (from a cluster run for example) and runs pipeline steps 1 and 2 on multiple threads
3) `update_github.sh` name says it all