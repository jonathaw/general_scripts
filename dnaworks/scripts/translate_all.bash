#!/bin/bash

input=$1
count=1;

pwdir=`pwd`; 
for s in `awk '{print $2}' $input `; do 
	./make_dnaworksinput.bash $s > job; 
	$pwdir/../bin/dnaworks.linux.f95release job ; 
	mv dnaworks.out $input.$count.dnaworks.out ; 
	./getSeq.pl $input.$count.dnaworks.out > seq; 
	mv seq $input.$count.seq;  
	rm job
	let count=count+1;
	#let $count=$count+1;  
done
