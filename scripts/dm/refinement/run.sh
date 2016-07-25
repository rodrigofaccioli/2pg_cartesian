#!/bin/bash


execution_list=$1
total_rep=$2
gmx_path=$3

min_mdp_all="min_all_implicit.mdp"
min_mdp_none="min_none_implicit.mdp"
eq_mdp="eq_implicit.mdp"
md_mdp="md_implicit.mdp"

total=$(wc -l "$execution_list" | awk '{print $1}')

n=1
while [ $n -le $total ]; do

	name=$(head -n $n "$execution_list" | tail -n 1 | sed 's/.pdb//g')

	if [ ! -e "$name"".md.1.gro" ]; then	# if this MD has not been made before...
		error=0

		"$gmx_path""./"pdb2gmx -f "$name"".pdb" -o "$name"".gro" -p "$name"".top" -ignh -water "none" -ff "amber99sb-ildn"
		if [ ! -e "$name"".gro" ]; then
			error=1
		fi

		"$gmx_path""./"editconf -f "$name"".gro" -o "$name"".gro" -c

		"$gmx_path""./"grompp -f "$min_mdp_none" -c "$name"".gro" -p "$name"".top" -o "$name"".min_none.tpr"

        "$gmx_path""./"mdrun -s "$name"".min_none.tpr" -v -deffnm "$name"".min_none" -pd

		"$gmx_path""./"grompp -f "$min_mdp_all" -c "$name"".min_none.gro" -p "$name"".top" -o "$name"".min_all.tpr"

		"$gmx_path""./"mdrun -s "$name"".min_all.tpr" -v -deffnm "$name"".min_all" -pd

		rep=1
		while [ $rep -le $total_rep ]; do

			"$gmx_path""./"grompp -f "$eq_mdp" -c "$name"".min_all.gro" -p "$name"".top" -o "$name"".eq.""$rep"".tpr"

                	"$gmx_path""./"mdrun -s "$name"".eq.""$rep"".tpr" -v -deffnm "$name"".eq.""$rep" -pd
	                if [ ! -e "$name"".eq.""$rep"".gro" ]; then
        	                error=1
                	fi


			"$gmx_path""./"grompp -f "$md_mdp" -c "$name"".eq.""$rep"".gro" -p "$name"".top" -o "$name"".md.""$rep"".tpr"

			"$gmx_path""./"mdrun -s "$name"".md.""$rep"".tpr" -v -deffnm "$name"".md.""$rep" -pd
			if [ ! -e "$name"".md.""$rep"".gro" ]; then
				error=1
			fi

			echo "C-alpha Protein" | "$gmx_path""./"editconf -f "$name"".md.""$rep"".gro" -o "$name"".md.""$rep""_last.pdb"

			echo "C-alpha Protein" | "$gmx_path""./"trjconv -f "$name"".md.""$rep"".xtc" -s "$name"".md.""$rep"".tpr" -o "$name"".md.""$rep"".pdb" -fit "rot+trans"

			echo "C-alpha C-alpha" | "$gmx_path""./"g_rms -f "$name"".md.""$rep"".xtc" -s "$name"".md.""$rep"".tpr" -o "$name"".rmsd_ca-ca_ini.""$rep"".xvg"

                	echo "Protein Protein" | "$gmx_path""./"g_hbond -f "$name"".md.""$rep"".xtc" -s "$name"".md.""$rep"".tpr" -num "$name"".hbond_prot-prot.""$rep"".xvg"

	                echo "MainChain+H MainChain+H" | "$gmx_path""./"g_hbond -f "$name"".md.""$rep"".xtc" -s "$name"".md.""$rep"".tpr" -num "$name"".hbond_main-main.""$rep"".xvg"

        	        echo "Protein Protein" | "$gmx_path""./"g_sas -f "$name"".md.""$rep"".xtc" -s "$name"".md.""$rep"".tpr" -o "temporary_sasa.xvg"
			grep "@" "temporary_sasa.xvg" | grep -v "@ s" > "$name"".aSASA.""$rep"".xvg"
			grep "@" "temporary_sasa.xvg" | grep -v "@ s" > "$name"".tSASA.""$rep"".xvg"
			grep "@" "temporary_sasa.xvg" | grep -v "@ s" > "$name"".tSASA.""$rep"".xvg"
			grep "@" "temporary_sasa.xvg" | grep -v "@ s" > "$name"".dG_solv.""$rep"".xvg"
			grep -v "#" "temporary_sasa.xvg" | grep -v "@" | awk '{print $1,$2}' >> "$name"".aSASA.""$rep"".xvg"
			grep -v "#" "temporary_sasa.xvg" | grep -v "@" | awk '{print $1,$3}' >> "$name"".pSASA.""$rep"".xvg"
			grep -v "#" "temporary_sasa.xvg" | grep -v "@" | awk '{print $1,$4}' >> "$name"".tSASA.""$rep"".xvg"
			grep -v "#" "temporary_sasa.xvg" | grep -v "@" | awk '{print $1,$5}' >> "$name"".dG_solv.""$rep"".xvg"
			rm "temporary_sasa.xvg"

                	echo "C-alpha C-alpha" | "$gmx_path""./"g_gyrate -f "$name"".md.""$rep"".xtc" -s "$name"".md.""$rep"".tpr" -o "$name"".rg-ca.""$rep"".xvg"

	                echo "Potential" | "$gmx_path""./"g_energy -f "$name"".md.""$rep"".edr" -o "$name"".potential.""$rep"".xvg"
        	        echo "GB-Polarization Nonpolar-Sol." | "$gmx_path""./"g_energy -f "$name"".md.""$rep"".edr" -sum -o "$name"".gbsa.""$rep"".xvg"
                	echo "LJ-14 LJ-(SR)" | "$gmx_path""./"g_energy -f "$name"".md.""$rep"".edr" -sum -o "$name"".vdw.""$rep"".xvg"
	                echo "Coulomb-14 Coulomb-(SR)" | "$gmx_path""./"g_energy -f "$name"".md.""$rep"".edr" -sum -o "$name"".coulomb.""$rep"".xvg"

			let rep=$rep+1

			rm \#*
		done

		if [ $error -eq 1 ]; then
			mkdir error/ 2>/dev/null
			cp "$name"".pdb" error/
		fi

	fi

	let n=$n+1
done
