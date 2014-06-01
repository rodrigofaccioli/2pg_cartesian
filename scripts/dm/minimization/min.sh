#!/bin/bash

# Executes minimization process of energy in all pdb files that are in current directory
#
# ./min <GROMACS path> 

gmx_path=$1

min_mdp_none="energy_minimization_implicit.mdp"
export GMX_MAXBACKUP=-1 # nao criar arquivos de backup do gromacs

here=$(pwd)
total_rep=1


#Identify how many pdbs there is in directory "$here"
cd "$here"
find . -maxdepth 1 -type f -name "*.pdb" | sed 's/.\///g' > "$here""/""temporary_list_pdbs"
total_pdbs=$(wc -l "temporary_list_pdbs" | awk '{print $1}')



#How many execution?
n=1
while [ $n -le $total_pdbs ]; do

	#get name of file and remove .pdb
	name=$(head -n $n "temporary_list_pdbs" | tail -n 1 | sed 's/.pdb//g')

	"$gmx_path""./"pdb2gmx -f "$name".pdb -o min.gro -p topol.top -water none -ff amber99sb-ildn -ignh
	"$gmx_path""./"grompp -f energy_minimization_implicit.mdp -c min.gro -p topol.top -o min_$n.tpr
	"$gmx_path""./"mdrun -s min_$n.tpr -v -deffnm min_$n -pd 
	"$gmx_path""./"editconf -f min_$n.gro -o "$name".pdb

	rm min_$n.tpr
	rm min_$n.gro
	rm min_$n.trr
	rm min_$n.log
	rm min_$n.edr

	let n=$n+1
done

rm mdout.mdp
rm min.gro
rm posre.itp
rm topol.top
rm temporary_list_pdbs

