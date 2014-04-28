#!/bin/bash

#----------------------------------------------------------------------------------
#
#     Script para rodar predições em cluster
#
#
#
#                                                         Leandro Oliveira Bortot
#                                                                       13 Sep 2013
#--------------------------------------------------------------------------------


script_evolutivo=$1
arq_parametros=$2
nodelist=$3
local_node=$4	# pasta nos nós...


pasta=$(pwd)
local="$pasta""/"



total_nodes=$(cat "$nodelist" | grep ";" | wc -l | awk '{print $1}')	# numero total de nós
let total_nodes=$total_nodes-1

node=1
while [ $node -le $total_nodes ]; do	# para cada nó...



	let linha=$node+1	# +1 devido ao cabeçalho no arquivo nodelist

	node_name=$(head -n $linha "$nodelist" | tail -n 1 | awk '{print $1}')
	total_cores=$(head -n $linha "$nodelist" | tail -n 1 | awk '{print $2}')
	echo -ne "$node_name... "

	cat "$arq_parametros" | grep "$node_name" | awk '{$1="" ; print $0}' | cut -c 2-9999  >> "temporario_""$arq_parametros""_""$node_name"

	total_pred_node=$(cat "temporario_""$arq_parametros""_""$node_name" | grep ";" | wc -l | awk '{print $1}')	# total de predições para o nó atual


	# Montagem dos arquivos de parametros para cada core individual
	i=1
	core=1
	while [ $i -le $total_pred_node ]; do	# para cada predição do nó...
		if [ $core -gt $total_cores ]; then
			core=1	# nao permite que o numero total de cores do nó seja ultrapassado
		fi
		
		head -n $i "temporario_""$arq_parametros""_""$node_name" | tail -n 1 >> "temporario_""$arq_parametros""_""$node_name""_""$core"

		let core=$core+1
		let i=$i+1
	done


	scp "$script_evolutivo" "$node_name":"$local_node" >/dev/null
	scp -r "arquivos/" "$node_name":"$local_node" >/dev/null

	total_submissoes_no=$(ls -lh "temporario_""$arq_parametros""_""$node_name""_"* | wc -l | awk '{print $1}')	# Total de submissoes por nó

	sub=1
	while [ $sub -le $total_submissoes_no ]; do

		head -n 1 "$arq_parametros" | awk '{$1="" ; print $0}' | cut -c 2-9999 > "$arq_parametros""_""$node_name""_""$sub"
		cat "temporario_""$arq_parametros""_""$node_name""_""$sub" >> "$arq_parametros""_""$node_name""_""$sub"


		scp "$arq_parametros""_""$node_name""_""$sub" "$node_name":"$local_node" >/dev/null
		nohup ssh "$node_name" "cd ""$local_node"" ; ./""$script_evolutivo"" ""$arq_parametros""_""$node_name""_""$sub""; exit "  > /dev/null 2>&1&
		cd "$local"

		let sub=$sub+1
	done


	echo "OK"


	let node=$node+1
done

rm temporario_*
