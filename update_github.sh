#!/bin/bash

file1="./mainscript.pl"
file2="./lib/*"
file3="README.md"

git commit -m 'auto commit' $file1
git commit -m 'auto commit' $file2
git commit -m 'auto commit' $file3

git push origin main
