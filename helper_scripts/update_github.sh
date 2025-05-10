#!/bin/bash

file1="../mainscript.pl"
file4="../README.md"
file5="./update_github.sh"
file6="../images/*"

git commit -m 'auto commit' $file1
git commit -m 'auto commit' $file4
git commit -m 'auto commit' $file5
git commit -m 'auto commit' $file6

git push origin main
