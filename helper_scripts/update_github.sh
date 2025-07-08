#!/bin/bash

file1="/home/peter/TE-discovery/mainscript.pl"
file4="/home/peter/TE-discovery/README.md"
file5="/home/peter/TE-discovery/helper_scritps/update_github.sh"
file6="/home/peter/TE-discovery/images/*"

git commit -m 'auto commit' $file1
git commit -m 'auto commit' $file4
git commit -m 'auto commit' $file5
git commit -m 'auto commit' $file6

git push origin main
