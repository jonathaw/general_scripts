#!/bin/bash

file=$1
#lines=`wc j| awk '{print $1-1 }'`

fiveprime="" #eg. CCTGATAAAGAAGCTTACTCAGGTTTT" 
threeprime="" #eg. GATGCATTGAAGCATGGTTTCGAA

rm seq
rm tmp
rm j
rm tmp2

echo $fiveprime > tmp
echo $threeprime > tmp2
awk '$2=="DNA"{print NR}' $file > j
solution=`head -1 j`

awk '{printf $1}' tmp > seq 
awk '{if (NR == '$solution' + 2) printf $2 }' $file >> seq
awk '{if (NR == '$solution' + 3) printf $2 }' $file >> seq
awk '{printf $1}' tmp2 >> seq
echo >> seq

#rm tmp j tmp2
#${var##[:space:]*
rm tmp tmp2
