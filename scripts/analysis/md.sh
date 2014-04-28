#!/bin/bash

#
#
#		Script para fazer dinamica molecular dos individuos não-dominados da última fronteira
#
#							Leandro Oliveira Bortot
#							21 Sep 2013
#
#

protpred=$1
gmx=$2
total_ind=$3
ger=$4
arq_config=$5
nt=$6
nativa=$7


dom=0	# Número de indivíduos que dominam os indivíduos da fronteira de intersse. dom=0 dá os indivíduos não-dominados, mas isso pode ser mudado.
avg_time_i=80	# tempo inicial para tomar a média da energia e RMSD
avg_time_f=100	# tempo final para tomar a média da energia e RMSD


saida_minimizacao_potencial="plot_min_pot.xvg"
saida_dinamica_potencial="plot_md_pot.xvg"
saida_dinamica_rmsd_ini="plot_md_rmsd-ini.xvg"
saida_dinamica_rmsd_nat="plot_md_rmsd-nat.xvg"

rm "$saida_minimizacao_potencial" 2>/dev/null
rm "$saida_dinamica_potencial" 2>/dev/null
rm "$saida_dinamica_rmsd_ini" 2>/dev/null
rm "$saida_dinamica_rmsd_nat" 2>/dev/null
rm "PROT_IND_FINAL_"*".pdb" 2>/dev/null
rm "temporario"* 2>/dev/null
rm -r "md/" 2>/dev/null
rm "front_""$ger""_""$dom""_"*".pdb" 2>/dev/null





echo -ne "\tGerando PDBs dos individuos da ultima geracao... "

"$protpred""src/"./protpred-Gromacs_Nerf_Population "$arq_config"  "pop_""$ger"".pop" >/dev/null 2>/dev/null	# Gera os PDBs dos indivíduos da geração-alvo

echo "OK"



echo -ne "\tPadronizando estrutura nativa..."
"$gmx"./pdb2gmx -f "$nativa" -o "nativa.pdb" -ignh -ff "charmm27" -water "none" >/dev/null 2>/dev/null
"$gmx"./editconf -f "nativa.pdb" -o "nativa.pdb" -resnr 1 >/dev/null 2>/dev/null
echo "OK"




echo -ne "\tSelecionando individuos não-dominados... "

cat "front_""$ger"".front" | grep -v ";" | awk '{print $3}' > "temporario_is_dominated_by"
total_front=0
ind=0	# numeração do 2pg começa com 0
while [ $ind -lt $total_ind ]; do	# para cada indivíduo ...

	let linha=$ind+1

	is_dominated_by=$(head -n $linha "temporario_is_dominated_by" | tail -n 1)

	if [ $is_dominated_by -eq $dom ]; then	# se o indivíduo é dominado pelo mesmo número dos de interesse (dom)

		cp "PROT_IND_FINAL_""$ind"".pdb" "front_""$ger""_""$dom""_""$ind"".pdb"	# copia o pdb do ind

		let total_front=$total_front+1	# incrementa o contador de ind não dominados na fronteira
	fi

	let ind=$ind+1
done

rm "temporario_is_dominated_by" "PROT_IND_FINAL_"*
echo "OK"



echo -e "\t""$total_front"" individuos nao-dominados foram encontrados na geracao ""$ger""."



mkdir md >/dev/null 2>/dev/null
cd md
mv ../"front_""$ger""_""$dom""_"*".pdb" .
cp ../*mdp .

	ind=1
	while [ $ind -le $total_front ]; do	# Para cada ind não-dominado...

		echo -e "\t\tSimulacao de dinamica molecular do individuo $ind de $total_front"


		"$gmx""/"./pdb2gmx -f "front_""$ger""_""$dom""_""$ind"".pdb" -o "$ind"".gro" -p "$ind"".top" -ignh -ff "charmm27" -water "none" >/dev/null 2>/dev/null

		echo -ne "\t\t\tMinimizacao de energia... "
		"$gmx""/"./grompp -f "min_none.mdp" -c "$ind"".gro" -p "$ind"".top" -o "$ind""_min_none.tpr" >/dev/null 2>/dev/null
		"$gmx""/"./mdrun -s "$ind""_min_none.tpr" -deffnm "$ind""_min_none" -nt $nt >/dev/null 2>/dev/null

		"$gmx""/"./grompp -f "min_all.mdp" -c "$ind""_min_none.gro" -p "$ind"".top" -o "$ind""_min_all.tpr" >/dev/null 2>/dev/null
		"$gmx""/"./mdrun -s "$ind""_min_all.tpr" -deffnm "$ind""_min_all" -nt $nt >/dev/null 2>/dev/null

		if [ -e "$ind""_min_all.gro" ]; then
			echo "potential" | "$gmx""/"./g_energy -f "$ind""_min_all.edr" -o "$ind""_min_all.xvg" >/dev/null 2>/dev/null
			pot=$(tail -n 1 "$ind""_min_all.xvg" | awk '{print $2}')
			echo "$ind"" ""$pot" >> "$saida_minimizacao_potencial"
			echo "OK"
		else
			echo "$ind"" " >> "$saida_minimizacao_potencial"
			echo "Erro"
		fi




		echo -ne "\t\t\tTermalizacao... "
		"$gmx""/"./grompp -f "eq.mdp" -c "$ind""_min_all.gro" -p "$ind"".top" -o "$ind""_eq.tpr" >/dev/null 2>/dev/null
		"$gmx""/"./mdrun -s "$ind""_eq.tpr" -deffnm "$ind""_eq" -nt $nt -v >/dev/null 2>/dev/null
		if [ -e "$ind""_eq.gro" ]; then
			echo "OK"
		else
			echo "Erro"
		fi


		echo -ne "\t\t\tDinamica... "
		"$gmx""/"./grompp -f "md.mdp" -c "$ind""_eq.gro" -p "$ind"".top" -o "$ind""_md.tpr" >/dev/null 2>/dev/null
		"$gmx""/"./mdrun -s "$ind""_md.tpr" -deffnm "$ind""_md" -nt $nt >/dev/null 2>/dev/null
		if [ -e "$ind""_eq.gro" ]; then


			echo "potential" | "$gmx""/"./g_energy -f "$ind""_md.edr" -o "$ind""_md_pot.xvg" -b $avg_time_i -e $avg_time_f >/dev/null 2>/dev/null
			cat "$ind""_md_pot.xvg" | grep -v "#" | grep -v "@" | grep -v ";" | awk '{print $2}' > "temporario_pot"
			total=$(wc -l "temporario_pot" | awk '{print $1}')			
			pot=$("$protpred""/scripts/analysis/"./avg-sd "temporario_pot" $total)

			echo "C-alpha C-alpha" | "$gmx""/"./g_rms -f "$ind""_md.xtc" -s "$ind""_md.tpr" -o "$ind""_md_rmsd-ini.xvg" -b $avg_time_i -e $avg_time_f >/dev/null 2>/dev/null
			cat "$ind""_md_rmsd-ini.xvg" | grep -v "#" | grep -v "@" | grep -v ";" | awk '{print $2}' > "temporario_rmsd-ini"
			total=$(wc -l "temporario_rmsd-ini" | awk '{print $1}')			
			rmsd_ini=$("$protpred""/scripts/analysis/"./avg-sd "temporario_rmsd-ini" $total)

			echo "C-alpha C-alpha" | "$gmx""/"./g_rms -f "$ind""_md.xtc" -s "../nativa.pdb" -o "$ind""_md_rmsd-nat.xvg" -b $avg_time_i -e $avg_time_f >/dev/null 2>/dev/null
			cat "$ind""_md_rmsd-nat.xvg" | grep -v "#" | grep -v "@" | grep -v ";" | awk '{print $2}' > "temporario_rmsd-nat"
			total=$(wc -l "temporario_rmsd-nat" | awk '{print $1}')			
			rmsd_nat=$("$protpred""/scripts/analysis/"./avg-sd "temporario_rmsd-nat" $total)

			echo "OK"

			echo "$ind"" ""$pot" >> "$saida_dinamica_potencial"
			echo "$ind"" ""$rmsd_ini" >> "$saida_dinamica_rmsd_ini"
			echo "$ind"" ""$rmsd_nat" >> "$saida_dinamica_rmsd_nat"
		else
			echo "Erro"

			echo "$ind"" " >> "$saida_dinamica_potencial"
			echo "$ind"" " >> "$saida_dinamica_rmsd_ini"
			echo "$ind"" " >> "$saida_dinamica_rmsd_nat"
		fi

		rm \#* temporario* 2>/dev/null

		let ind=$ind+1
	done

mv plot_*xvg ../
cd ../

