#!bin/bash
for file in $1/*.bed; do
    bedtools getfasta -fi ../dm3.fa -bed $file -fo $file.fa
done
