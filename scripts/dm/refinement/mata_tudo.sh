#!/bin/bash

nodefile=$1
path_to_pdbs=$2

here=$(pwd)
total_rep=5


total_nodes=$(wc -l "$nodefile" | awk '{print $1}')

node=1
while [ $node -le $total_nodes ]; do

	node_name=$(head -n $node "$nodefile" | tail -n 1)

	echo -n "$node_name"" ... "

	ssh "$node_name" "killall run.sh"
	ssh "$node_name" "killall mdrun"
	ssh "$node_name" "killall run.sh"

	ssh "$node_name" "rm -r ""$path_to_pdbs""/"

	echo "ok"

	let node=$node+1
done
