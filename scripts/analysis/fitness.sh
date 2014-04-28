#!/bin/bash

#----------------------------------------------------------------------------------
#
#     Script para calcular o comportamento dos valores de fitness dos objetivos ao longo das gerações
#
#	
#
#                                                                       25 jan 2012
#--------------------------------------------------------------------------------

arq_fitness=$1
total_ind=$2
total_ger=$3

a=0
b=0
linhas=0

# titulo será inserido nos nomes dos arquivos para diferenciar os fitness
titulo=$(echo $arq_fitness | sed 's/.fit//g')

# se arquivos de saída existem, são apagados
if [ -e fitness_"$titulo"_fitness ]; then rm fitness_"$titulo"_fitness; fi
if [ -e fitness_"$titulo"_ger ]; then rm fitness_"$titulo"_ger; fi
if [ -e fitness_"$titulo"_maior ]; then rm fitness_"$titulo"_maior; fi
if [ -e fitness_"$titulo"_menor ]; then rm fitness_"$titulo"_menor; fi
if [ -e fitness_"$titulo"_maior_log ]; then rm fitness_"$titulo"_maior_log; fi
if [ -e fitness_"$titulo"_menor_log ]; then rm fitness_"$titulo"_menor_log; fi

# cria os arquivos de saída vazios
touch fitness_"$titulo"_menor fitness_"$titulo"_maior fitness_"$titulo"_maior_log fitness_"$titulo"_menor_log temp_ger temp_fitness

i=1

while [ $i -le $total_ger ]; do

	# Definição da linha a ser lida no arquivo de fitness. a-> contribuição dos indivíduos; b-> contribuição do cabelhaço de cada geração
	let a=$total_ind*$i
	let b=2*$i
	let linhas=$a+$b

	# Separa os valores de fitness da geração atual em um arquivo temporário
	head -n $linhas $arq_fitness | tail -n $total_ind | sed 's/\ /\ \ \ \ \ \ \ \ \ /g' | cut -c 10-100 | sed 's/\ //g' > temporario_fitness_"$titulo"

	
	j=2
	while [ $j -le $total_ind ]; do

		# menor e maior sao, ambos, o primeiro valor
		menor=$(head -n 1 temporario_fitness_"$titulo")
		maior=$(head -n 1 temporario_fitness_"$titulo")

		# valor começa na 2 linha (j=2)
		valor=$(head -n $j temporario_fitness_"$titulo" | tail -n 1)

		# faz comparaçoes. 1=verdadeiro, 0=falso
		compara_menor=$(echo "$valor"" < ""$menor" | bc)
		compara_maior=$(echo "$valor"" > ""$maior" | bc)

		# troca valores se necessário
		if [ $compara_menor -eq 1 ]; then menor=$valor; fi
		if [ $compara_maior -eq 1 ]; then maior=$valor; fi

		let j=$j+1

	done

	maior_log=$(echo "l(""$maior"")" | bc -l)
	maior_log=$(echo "$maior_log"" / 2.302585093" | bc -l)

        menor_log=$(echo "l(""$menor"")" | bc -l)
        menor_log=$(echo "$menor_log"" / 2.302585093" | bc -l)


	echo $i $maior >> fitness_"$titulo"_maior
	echo $i $menor >> fitness_"$titulo"_menor

	echo $i $maior_log >> fitness_"$titulo"_maior_log
	echo $i $menor_log >> fitness_"$titulo"_menor_log

	echo $i

	let i=$i+1

done

i=1

while [ $i -le $total_ger ]; do
	j=1
	while [ $j -le $total_ind ]; do
		echo $i >> temp_ger
		let j=$j+1
	done
	let i=$i+1
done

paste temp_ger temp_fitness > grafico.xvg
