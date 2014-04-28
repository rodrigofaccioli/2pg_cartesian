#!/bin/bash

#
#
#	Script para analisar predições multi-objetivo. São analisados:
#		- A fronteira não-dominada de cada geração
#		- O valor de cada objetivo dos indivíduos não-dominados ao longo das gerações
#		- Número de indivíduos não-dominados por geração
#		- Valor de RMSD dos indivíduos não-dominados
#		- Valor de GDT-TS dos indivíduos não-dominados
#
#							Leandro Oliveira Bortot
#							21 Sep 2013
#
#

protpred=$1
gmx=$2
maxcluster=$3
total_ind=$4
total_ger=$5
total_obj=$6
nativa=$7
arq_config=$8
arq_obj=$9

dom=0	# Número de indivíduos que dominam os indivíduos da fronteira de intersse. dom=0 dá os indivíduos não-dominados, mas isso pode ser mudado.


# Apaga arquivos pre-existentes
if [ $(ls -lh "front_""$ger""_""$dom""_"*".pdb" 2>/dev/null | wc -l | awk '{print $1}') -ge 1 ]; then rm "front_""$ger""_""$dom""_"*".pdb" ; fi 
if [ $(ls -lh "PROT_IND_FINAL"*".pdb" 2>/dev/null | wc -l | awk '{print $1}') -ge 1 ]; then rm "PROT_IND_FINAL"*".pdb" ; fi
if [ $(ls -lh "temporario_"* 2>/dev/null | wc -l | awk '{print $1}') -ge 1 ]; then rm "temporario_"* ; fi  

if [ -e "front_""$ger""_""$dom"".xvg" ]; then rm "front_""$ger""_""$dom"".xvg" ; fi
if [ $(ls -lh "plot_"*"_front_""$dom""_scatter.xvg" 2>/dev/null | wc -l | awk '{print $1}') -ge 1 ]; then rm "plot_"*"_front_""$dom""_scatter.xvg" ; fi
if [ $(ls -lh "plot_"*"_front_""$dom""_avg-sd.xvg" 2>/dev/null | wc -l | awk '{print $1}') -ge 1 ]; then rm "plot_"*"_front_""$dom""_avg-sd.xvg" ; fi

if [ -e "plot_nfront_""$dom"".xvg" ]; then rm "plot_nfront_""$dom"".xvg" ; fi

if [ -e "plot_rmsd_scatter.xvg" ]; then rm "plot_rmsd_scatter.xvg" ; fi
if [ -e "plot_rmsd_min.xvg" ]; then rm "plot_rmsd_min.xvg" ; fi
if [ -e "plot_rmsd_avg-sd.xvg" ]; then rm "plot_rmsd_avg-sd.xvg" ; fi
if [ -e "rmsd_min.pdb" ]; then rm "rmsd_min.pdb" ; fi

if [ -e "plot_gdt-ts_scatter.xvg" ]; then rm "plot_gdt-ts_scatter.xvg" ; fi
if [ -e "plot_gdt-ts_max.xvg" ]; then rm "plot_gdt-ts_max.xvg" ; fi
if [ -e "plot_gdt-ts_avg-sd.xvg" ]; then rm "plot_gdt-ts_avg-sd.xvg" ; fi
if [ -e "gdt-ts_max.pdb" ]; then rm "gdt-ts_max.pdb" ; fi


echo -ne "\tPadronizando estrutura nativa..."
"$gmx"./pdb2gmx -f "$nativa" -o "nativa.gro" -ignh -ff "charmm27" -water "none" -renum >/dev/null 2>/dev/null
"$gmx"./pdb2gmx -f "$nativa" -o "nativa.pdb" -ignh -ff "charmm27" -water "none" -renum >/dev/null 2>/dev/null
echo "OK"


ger=1
while [ $ger -le $total_ger ]; do	# Para cada geração...

	echo -e "\tGeracao $ger"

	echo -ne "\t\tSelecionando individuos não-dominados... "
# Leitura dos valores dos objetivos de todos os indivíduos da geração
	junta=""	# string que será construída para o comando paste

	i=1
	while [ $i -le $total_obj ]; do	# para cada obj...
		obj=$(head -n $i "$arq_obj" | tail -n 1 | awk '{print $1}')	# lê os objetivos que estão no arq_obj

		cat "$obj""_""$ger"".fit" | grep -v ";" | awk '{print $2}' > "temporario_obj""$i"	# Lê os valores dos obj para a geração-alvo a partir do arquivo .fit correspondente

		junta="$junta"" temporario_obj""$i"	# aumenta a string que terá o comando paste

		let i=$i+1
	done

	paste $junta > "temporario_obj"	# esse arquivo temporário tem os valores de obj de todos os ind da geração-alvo


	
# Filtragem: Somente os indivíduos não-dominados são considerados	
	total_front=0
	i=1
	while [ $i -le $total_ind ]; do	# para cada indivíduo ...

		cat "front_""$ger"".front" | grep -v ";" | awk '{print $3}' > "temporario_is_dominated_by"	# Separa o número de ind que domina cada um
		is_dominated_by=$(head -n $i "temporario_is_dominated_by" | tail -n 1)	# lê o número de ind que dominam o atual

		if [ $dom -eq $is_dominated_by ]; then	# se o indivíduo é dominado pelo mesmo número dos de interesse (dom)
			head -n $i "temporario_obj" | tail -n 1 >> "temporario_front_obj"	# copia valores dos obj
			let ind=$i-1	# numeração dos ind no 2PG começa do 0
			echo $ind >> "temporario_front_ind"	# copia numero do ind
			let total_front=$total_front+1	# incrementa o contador de ind não dominados na fronteira
		fi

		let i=$i+1
	done

	echo "OK"
	
	echo -ne "\t\tSelecionando os valores de objetivo para os individuos nao-dominados... "
	
# Ordenação pelo valor do primeiro objetivo
	paste "temporario_front_obj" "temporario_front_ind" > "temporario_front"	# junta objetivos e numero dos ind para a ordenação
	sort -g "temporario_front" > "temporario_front_sorted"	# ordena pela 1 coluna (primeiro objetivo)
	
	let colunas=$total_obj+1	# numero de colunas do arquivo temporario_front (+1 pq do numero do ind)
	cat "temporario_front_sorted" |  awk '{print $'"$colunas"'}' > "temporario_ind_sorted"	# awk dá a coluna $colunas, que é a última e portanto a que tem o número de cada ind

	
	
# Plotagem da fronteira
	cat "temporario_front_sorted" |  awk '{out=$1 ; for(i=2;i<='"$total_obj"';i++) out=out" "$i; print out}' > "front_""$ger""_""$dom"".xvg"	# awk dá as colunas de 1 até total_obj, ou seja, os valores dos objetivos crescentemente ordenados pelo primeiro


# Plota os valores dos objetivos para os indivíduos da fronteira
	i=1
	while [ $i -le $total_front ]; do	# para cada ind da fronteira...
		echo "$ger" >> "temporario_ger"	# cria arquivo temporário com o número da geração para cada ind da fronteira. Será usado para plotar os valores dos objetivos no próximo passo		
		let i=$i+1
	done
	
	
	
	i=1
	while [ $i -le $total_obj ]; do

		obj=$(head -n $i "$arq_obj" | tail -n 1 | awk '{print $1}')	# lê os objetivos que estão no arq_obj

		cat "temporario_front_sorted" |  awk '{print $'"$i"'}' > "temporario_""$obj""_front_""$ger""_""$dom"".xvg"	# awk dá as colunas de 1 até total_obj, ou seja, os valores dos objetivos crescentemente ordenados pelo primeiro

		paste "temporario_ger" "temporario_""$obj""_front_""$ger""_""$dom"".xvg" >> "plot_""$obj""_front_""$dom""_scatter.xvg"

		avgsd=$("$protpred""scripts/analysis/"./avg-sd "temporario_""$obj""_front_""$ger""_""$dom"".xvg" $total_front)
		echo -e "$ger""\t""$avgsd" >> "plot_""$obj""_front_""$dom""_avg-sd.xvg"

		let i=$i+1
	done




# Plota o número de ind na fronteira
	echo "$ger $total_front" >> "plot_nfront_""$dom"".xvg"



	echo "OK"




	echo -ne "\t\tSalvando PDBs dos individuos nao-dominados... "


# PDBs dos indivíduos não-dominados
	"$protpred""src/"./protpred-Gromacs_Nerf_Population "$arq_config"  "pop_""$ger"".pop" >/dev/null 2>/dev/null	# Gera os PDBs dos indivíduos da geração-alvo

	i=1
	while [ $i -le $total_front ]; do	# para cada ind da fronteira...
		ind=$(head -n $i "temporario_ind_sorted" | tail -n 1)	# lê número do ind
		cp "PROT_IND_FINAL_""$ind"".pdb" "front_""$ger""_""$dom""_""$i""_""$ind"".pdb"	# copia o pdb do ind
		
		echo "$ger" >> "temporario_ger"	# cria arquivo temporário com o número da geração para cada ind da fronteira. Será usado para plotar os valores dos objetivos no próximo passo
		
		let i=$i+1
	done
	rm "PROT_IND_FINAL_"*".pdb"
	
	echo "OK"
	
	echo -ne "\t\tCalculando RMSDs e selecionando o menor... "

# RMSD
	i=1
	while [ $i -le $total_front ]; do	# Para cada ind não-dominado...

		"$gmx"./pdb2gmx -f "front_""$ger""_0_""$i""_"* -o "temporario.gro" -ignh -ff "charmm27" -water "none" -renum >/dev/null 2>/dev/null
		echo "C-alpha C-alpha" | "$gmx"./g_rms -f "temporario.gro" -s "nativa.gro" -o "rmsd.xvg" >/dev/null 2>/dev/null	# Calcula RMSD em relação à nativa

		valor=$(tail -n 1 "rmsd.xvg" | awk '{print $2}')	# Extrai valor do RMSD

		echo -e "$ger""\t""$valor" >> "plot_rmsd_scatter.xvg"	# Plot scatter com os RMSDs

		echo -e "$valor""\t""$i" >> "temporario_rmsd"	# Temporário para encontrar o RMSD minimo
		echo "$valor" >> "temporario_rmsd_avg-sd"

		rm "rmsd.xvg" "temporario.gro" \#* 2>/dev/null

		let i=$i+1

	done


	avgsd=$("$protpred""scripts/analysis/"./avg-sd "temporario_rmsd_avg-sd" $total_front)
	echo -e "$ger""\t""$avgsd" >> "plot_rmsd_avg-sd.xvg"

	sort -n "temporario_rmsd" > "temporario_rmsd_sorted"	# Ordena para encontrar o RMSD mínimo
	rmsd=$(head -n 1 "temporario_rmsd_sorted" | awk '{print $1}')	# Valor do menor RMSD
	ind=$(head -n 1 "temporario_rmsd_sorted" | awk '{print $2}')	# Indivíduo que tem o menor RMSD
	echo -e "$ger""\t""$rmsd" >> "plot_rmsd_min.xvg"	# Plot de RMSD mínimo

	# Trajetória com os indivíduos de menor RMSD
	cp "front_""$ger""_0_""$ind""_"* "temporario.pdb"
	"$gmx"./pdb2gmx -f "temporario.pdb" -o "temporario.pdb" -ignh -ff "charmm27" -water "none" -renum >/dev/null 2>/dev/null
	echo "MODEL ""$ger" >> "rmsd_min.pdb"
	cat "temporario.pdb" | grep "ATOM" >> "rmsd_min.pdb"
	echo "ENDMDL" >> "rmsd_min.pdb"
	rm "temporario.pdb" \#* 2>/dev/null

	
	echo "OK"
	echo -ne "\t\tCalculando GDT-TSs e selecionando o maior... "
	
	
# GDT-TS
	i=1
	while [ $i -le $total_front ]; do	# Para cada ind não-dominado...
		"$gmx"./pdb2gmx -f "front_""$ger""_0_""$i""_"* -o "temporario.pdb" -ignh -ff "charmm27" -water "none" -renum >/dev/null 2>/dev/null

		valor=$("$maxcluster"./maxcluster -e "nativa.pdb" -p "temporario.pdb" -gdt | tail -n 1 | awk '{print $2}')

		echo -e "$ger""\t""$valor" >> "plot_gdt-ts_scatter.xvg"	# Plot scatter com os GDT-TS

		echo -e "$valor""\t""$i" >> "temporario_gdt-ts"	# Temporário para encontrar o GDT-TS máximo
		echo "$valor" >> "temporario_gdt-ts_avg-sd"

		rm  \#* 2>/dev/null
		
		let i=$i+1
	done

	
	avgsd=$("$protpred""scripts/analysis/"./avg-sd "temporario_gdt-ts_avg-sd" $total_front)
	echo -e "$ger""\t""$avgsd" >> "plot_gdt-ts_avg-sd.xvg"
	
	sort -nr "temporario_gdt-ts" > "temporario_gdt-ts_sorted"	# Ordena para encontrar o GDT-TS máximo
	gdtts=$(head -n 1 "temporario_gdt-ts_sorted" | awk '{print $1}')	# Valor do maior GDT-TS
	ind=$(head -n 1 "temporario_gdt-ts_sorted" | awk '{print $2}')	# Indivíduo que tem o maior GDT-TS
	echo -e "$ger""\t""$gdtts" >> "plot_gdt-ts_max.xvg"	# Plot de GDT-TS máximo

	# Trajetória com os indivíduos de maior GDT-TS
	cp "front_""$ger""_0_""$ind""_"* "temporario.pdb"
	"$gmx"./pdb2gmx -f "temporario.pdb" -o "temporario.pdb" -ignh -ff "charmm27" -water "none" -renum >/dev/null 2>/dev/null
	echo "MODEL ""$ger" >> "gdt-ts_max.pdb"
	cat "temporario.pdb" | grep "ATOM" >> "gdt-ts_max.pdb"
	echo "ENDMDL" >> "gdt-ts_max.pdb"
	rm "temporario.pdb" \#* 2>/dev/null
	
	
	echo "OK"
	

	let ger=$ger+1
	rm temporario_*
	rm "front_"*"_""$dom""_"*"_"*".pdb"

done

echo
echo "\tFim das geracoes"
echo

echo -ne "\tAlinhando trajetoria dos melhores RMSD e GDT-TS... "

mv "rmsd_min.pdb" "temporario_rmsd_min.pdb"
echo "C-alpha Protein" | "$gmx"./trjconv -f "temporario_rmsd_min.pdb" -s "nativa.gro" -fit "rot+trans" -o "rmsd_min.pdb" >/dev/null 2>/dev/null

mv "gdt-ts_max.pdb" "temporario_gdt-ts_max.pdb"
echo "C-alpha Protein" | "$gmx"./trjconv -f "temporario_gdt-ts_max.pdb" -s "nativa.gro" -fit "rot+trans" -o "gdt-ts_max.pdb" >/dev/null 2>/dev/null

echo "OK"

rm temporario_* \#*

