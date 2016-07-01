#!/bin/bash

rm tmp tmp2 myseq.fasta

for i in *seq; do echo ">$i" >> myseq.fasta; cat $i >> myseq.fasta  ;done 
mkdir out
gzip *out 
mv *out.gz out
mv *seq out

