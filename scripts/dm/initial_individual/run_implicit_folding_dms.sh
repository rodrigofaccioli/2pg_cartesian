#!/bin/bash


md_first=$1
md_last=$2
input_structure=$3	#optional if not doing initialization and minimization
native_structure=$4	#optional with analysis does not include native RMSD

gmx_path="/usr/local/bin"
mdp_minimization="min_all_implicit.mdp"
mdp_md="md_implicit.mdp"

md_time="100"	#ps
timestep="0.002"	#ps
temperature="309"	#K

max_error=5	# max times a MD can fail before trying to rerun

initial=0
minimization=0
mds=0
analysis=1

if [ $initial -eq 1 ]; then
	# Initial reading
	"$gmx_path"/./pdb2gmx -f "$input_structure" -o "protein.gro" -p "protein.top" -ignh -ff "amber99sb-ildn" -water "none"
fi

if [ $minimization -eq 1 ]; then
	# Minimization
	"$gmx_path"/./grompp -f "$mdp_minimization" -c "protein.gro" -p "protein.top" -o "minimization.tpr"

	"$gmx_path"/./mdrun -s "minimization.tpr" -v -deffnm "minimization" -nt 1

	rm \#*
fi

if [ $mds -eq 1 ]; then
	# Molecular dynamics

	# Determine the total number of 2fs steps and build a proper mdp file
	total_steps=$(echo "$md_time"" / ""$timestep" | bc)
	sed 's/_STEPS_/'"$total_steps"'/g' "$mdp_md" | sed 's/_TIMESTEP_/'"$timestep"'/g' | sed 's/_TEMPERATURE_/'"$temperature"'/g' > "temporary.mdp"

	n=$md_first
	error=0	# error counter

        while [ $n -le  $md_last ]; do
                # Verifies for MDs which were already made. This allow continuing executions
                if [ ! -e "md.""$n"".gro" ]; then

                        error=0
			# The executions are trapped inside an infinite loop. The only way to scape is via a break which happens if the DM was sucessfull or when it fails enough times.
                        while [ 1 -eq 1 ]; do
				"$gmx_path"/./grompp -f "temporary.mdp" -c "minimization.gro" -p "protein.top" -o "md.""$n"".tpr"
				"$gmx_path"/./mdrun -s "md.""$n"".tpr" -deffnm "md.""$n" -v -pd

				# If the output .gro exists, the DM was sucessfull. Continue to next one.
                                if [ -e "md.""$n"".gro" ]; then
                                        break
                                else
				# If the output .gro does not exist, the DM was not sucessfull. Increase error counter
                                        let error=$error+1      # increase the error counter

                                        # If this MD has already failed max_error times, give up
                                        if [ $error -ge $max_error ]; then
                                                break
                                        fi
                                fi
                        done
                fi
                let n=$n+1

		rm \#*
        done
fi


if [ $analysis -eq 1 ]; then

	rm "md.pdb" "temporary.pdb" >/dev/null 2>/dev/null

	n=$md_first
	while [ $n -le  $md_last ]; do

		if [ ! "$native_structure" == "" ]; then
			"$gmx_path"/./pdb2gmx -f "$native_structure" -o "native.gro" -p "native.top" -ignh -ff "amber99sb-ildn" -water "none"

			echo "C-alpha C-alpha" | "$gmx_path"/./g_rms -f "md.""$n"".xtc" -s "native.gro" -o "rmsd_ca-ca_nat.""$n"".xvg"
		fi

		echo "C-alpha C-alpha" | "$gmx_path"/./g_rms -f "md.""$n"".xtc" -s "md.""$n"".tpr" -o "rmsd_ca-ca_ini.""$n"".xvg"

		echo "Protein Protein" | "$gmx_path"/./g_hbond -f "md.""$n"".xtc" -s "md.""$n"".tpr" -num "hbond_prot-prot.""$n"".xvg"

		echo "MainChain+H MainChain+H" | "$gmx_path"/./g_hbond -f "md.""$n"".xtc" -s "md.""$n"".tpr" -num "hbond_main-main.""$n"".xvg"

		echo "Protein Protein" | "$gmx_path"/./g_sas -f "md.""$n"".xtc" -s "md.""$n"".tpr" -o "sasa.""$n"".xvg"
		grep -v "#" "sasa.""$n"".xvg" | grep -v "@" | awk '{print $1,$2}' > "aSASA.""$n"".xvg"
		grep -v "#" "sasa.""$n"".xvg" | grep -v "@" | awk '{print $1,$3}' > "pSASA.""$n"".xvg"
		grep -v "#" "sasa.""$n"".xvg" | grep -v "@" | awk '{print $1,$4}' > "tSASA.""$n"".xvg"
		grep -v "#" "sasa.""$n"".xvg" | grep -v "@" | awk '{print $1,$5}' > "dGsolv.""$n"".xvg"
		rm "sasa.""$n"".xvg"

		echo "C-alpha C-alpha" | "$gmx_path"/./g_gyrate -f "md.""$n"".xtc" -s "md.""$n"".tpr" -o "rg-ca.""$n"".xvg"

		echo "Potential" | "$gmx_path"/./g_energy -f "md.""$n"".edr" -o "potential.""$n"".xvg"
		echo -e "GB-Polarization""\n""Nonpolar-Sol." | "$gmx_path"/./g_energy -f "md.""$n"".edr" -sum -o "gbsa.""$n"".xvg"
		echo -e "LJ-14""\n""LJ-(SR)" | "$gmx_path"/./g_energy -f "md.""$n"".edr" -sum -o "vdw.""$n"".xvg"
		echo -e "Coulomb-14""\n""Coulomb-(SR)" | "$gmx_path"/./g_energy -f "md.""$n"".edr" -sum -o "coulomb.""$n"".xvg"

		echo "Protein" | "$gmx_path"/./editconf -f "md.""$n"".gro" -o "md.""$n"".pdb"
		echo "MODEL ""$n" >> "temporary.pdb"
		grep "ATOM" "md.""$n"".pdb" >> "temporary.pdb"
		echo "ENDMDL" >> "temporary.pdb"

		let n=$n+1

	done

	echo "C-alpha Protein" | "$gmx_path"/./trjconv -f "temporary.pdb" -s "temporary.pdb" -o "md.pdb" -fit "rot+trans"
	rm "temporary.pdb"

	rm \#*
fi
