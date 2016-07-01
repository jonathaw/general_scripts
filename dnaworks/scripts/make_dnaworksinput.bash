#!/bin/bash

seq=$1

echo "solutions 1"
echo "repeat 8"
echo "LOGFILE \"dnaworks.out\""

echo "pattern"
echo "BamHI GGATCC"
echo "XhoI  CTCGAG"
echo "NheI  GCTAGC"
echo "NdeI  CATATG"
echo "BsaI  GGTCTC"
echo "AarI CACCTGC"

echo "E. coli class II"

echo "protein"
echo $seq
echo "//"


