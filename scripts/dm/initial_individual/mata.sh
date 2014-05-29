#!/bin/bash

n=1

while [ $n -le 10 ]; do

	ssh vm"$n" "killall mdrun"
	ssh vm"$n" "killall run_implicit_folding_dms.sh"
	ssh vm"$n" "killall mdrun"
	ssh vm"$n" "killall run_implicit_folding_dms.sh"

	let n=$n+1

done



