#!/bin/bash

input_structure=$1
directory=$2

n=1
i=1
f=100
while [ $n -le 10 ]; do

	ssh vm"$n" "mkdir -p ""$directory""/"
	scp arquivos/* vm"$n":"$directory" >/dev/null 2>/dev/null
	nohup ssh vm"$n" "cd ""$directory"" ; ./run_implicit_folding_dms.sh $i $f ""$input_structure""" >/dev/null 2>&1&

	let n=$n+1
	let i=$i+100
	let f=$f+100

done



