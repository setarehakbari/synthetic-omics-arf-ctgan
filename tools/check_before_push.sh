#!/bin/bash
# Run from the repository root before every push:  bash tools/check_before_push.sh
# It lists (1) data-like files that would be committed, (2) local absolute paths, (3) file names with unusual characters.
echo "== data-like files tracked or untracked-but-not-ignored (should be empty except seeds.csv) =="
git ls-files --cached --others --exclude-standard | grep -E '\.(csv|tsv|xlsx|rds|RData|rda|pt|pkl|parquet|json)$'
echo "== local absolute paths in text files (should be empty) =="
grep -rIn --exclude-dir=.git --exclude-dir=tools -E '/Users/|/home/[a-z]|C:\\\\Users|Desktop' . | head -20
echo "== file names with characters other than letters, digits, . _ - / (should be empty) =="
git ls-files --cached --others --exclude-standard | grep -vE '^[A-Za-z0-9._/-]+$'
echo "done"
