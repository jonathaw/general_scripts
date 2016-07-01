#!/bin/bash

rm tmp tmp2 myseq.fasta
ls *seq  > j
n=`wc j | awk '{print $1}'`
echo $n
for i in `seq 1 $n`; do echo ">all.1line.$i.seq" >> myseq.fasta; cat all.1line.$i.seq >> myseq.fasta ;done 

#for i in *seq; do echo ">$i" >> myseq.fasta; cat $i >> myseq.fasta  ;done 
mkdir out
gzip *out 
mv *out.gz out
mv *seq out

