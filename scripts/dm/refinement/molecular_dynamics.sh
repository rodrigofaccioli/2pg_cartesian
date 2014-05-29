#!/bin/bash

nodefile=$1
path_to_pdbs=$2
gmxpath=$3

export GMX_MAXBACKUP=-1 # nao criar arquivos de backup do gromacs

here=$(pwd)
total_rep=1


cd "$path_to_pdbs"
find . -maxdepth 1 -type f -name "*.pdb" | sed 's/.\///g' > "$here""/""temporary_list_pdbs"
cd "$here"
total_pdbs=$(wc -l "temporary_list_pdbs" | awk '{print $1}')


total_nodes=$(wc -l "$nodefile" | awk '{print $1}')


#creates the execution list
rm "execution_list_"* 2>/dev/null
node=1
pdb=1
while [ $pdb -le $total_pdbs ]; do

	node_name=$(head -n $node "$nodefile" | tail -n 1)
	pdb_name=$(head -n $pdb "temporary_list_pdbs" | tail -n 1)

	echo "$pdb_name" >> "execution_list_""$node_name"

	if [ $node -eq $total_nodes ]; then
		node=1
	else
		let node=$node+1
	fi

	let pdb=$pdb+1
done


rm "temporary_list_pdbs"


node=1
while [ $node -le $total_nodes ]; do

	node_name=$(head -n $node "$nodefile" | tail -n 1)

	ssh "$node_name" mkdir -p "$path_to_pdbs"
	scp  "execution_list_""$node_name" "$node_name":"$path_to_pdbs""/"
	scp *".mdp" "$node_name":"$path_to_pdbs""/"
	scp *".sh" "$node_name":"$path_to_pdbs""/"
	scp "$path_to_pdbs""/"*".pdb" "$node_name":"$path_to_pdbs""/"
	nohup ssh "$node_name" "cd ""$path_to_pdbs"" ; ./run.sh execution_list_""$node_name"" ""$total_rep" "$gmxpath"  > /dev/null 2>&1&

	let node=$node+1
done
