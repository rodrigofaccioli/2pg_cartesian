#!/bin/bash

#----------------------------------------------------------------------------------
#
#     Script para rodar predições sequencialmente
#
#		Basta informar os parâmetros no arquivo parametros.config
#
#
#                                                         Leandro Oliveira Bortot
#                                                                       13 Sep 2013
#--------------------------------------------------------------------------------

arq_parametros=$1

pasta=$(pwd)
local="$pasta""/"

arq_config="configuracao.conf"
arquivos="$local"arquivos/	 #pasta que contém os FASTA, .pdb e .mdp

path_protpred="/home/leandro/programs/2pg/"
path_gromacs="/usr/local/gromacs/bin/"
path_maxcluster="/home/leandro/programs/"
objective_analisys_dimo_source="/home/faccioli/workspace/dimo/DIMO2"
objective_analisys="none" 


CRIA_ARQ_CONFIG=1 # criação do arquivo de configuração?

POP_INI=0 # criação da população inicial?

EA=0 # algoritmo evolutivo?

MONO_ANALYSIS=0		# Análise de algoritos Mono-Objetivo

FRONT_ANALYSIS=0	# Análise de algoritmo Multi-Objetivo

MOLECULAR_DYNAMICS=0	# Dinâmica Molecular nos indivíduos não-dominados da última geração

TAR=1	# Criar .tar.gz da pasta da predição no fim







# total_predicoes é o número de linhas do arquivo de parametros -1 devido ao cabeçalho. As linhas devem ter ";" no final, incluindo o cabeçalho
total_predicoes=$(cat $arq_parametros | grep ";" | wc -l |  awk '{print $1}')
let total_predicoes=$total_predicoes-1

# começa na predição 1
pred_atual=1

# enquanto pred_atual for menor ou igual a total_predicoes ...
while [ $pred_atual -le $total_predicoes ]; do

	# +1 para que a leitura dos parâmetros pule a primeira linha do arquivo, que tem o cabeçalho
	let linha=$pred_atual+1

	# leitura dos parâmetros...
	proteina=$(head -n $linha $arq_parametros | tail -n 1 | cut -c 1-4)
	titulo=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $1}')
	repeticoes=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $2}')
	N_inicial=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $3}')
	minimiza=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $4}')
	algoritmo=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $5}')
	nt=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $6}')
	objetivos=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $7}')
	geracoes=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $8}')
	individuos=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $9}')
	novosind=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $10}')
	archive=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $11}')
	crossover=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $12}')
	blx_alpha=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $13}')
	p_blx=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $14}')
	p_1pt=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $15}')
	p_2pt=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $16}')
	mutacao_ind=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $17}')
	mutacao=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $18}')
	max_mut=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $19}')

	par=20	# contem o número do próximo parâmetro a ser lido do arquivo
	i=1	# número do objetivo atual
	while [ $i -le $objetivos ]; do	# para cada objetivo
		obj[$i]=$(head -n $linha $arq_parametros | tail -n 1 | awk '{print $'"$par"'}')	# cada objetivo é armazenado em um elemento do vetor obj[]
		let i=$i+1
		let par=$par+1
	done


	rep=$N_inicial	# Repeticao atual é N_inicial, especificado no arquivo de parametros

	while [ $rep -le $repeticoes ]; do	# repete a mesma predicao
		
		# cria a pasta da simulacao
		mkdir "$titulo""_""$rep" 2>/dev/null

		# entra na pasta recem-criada
		cd "$titulo""_""$rep"

		# copia arquivos essenciais
		cp "$arquivos"* .

	# cria o arquivo de configuracao

		if [ $CRIA_ARQ_CONFIG -eq 1 ]; then

			if [ -e $arq_config ]; then rm $arq_config; fi

			echo "gromacs_energy_min = ""$minimiza" >> $arq_config
			echo "gromacs_energy_min_gen_oper = ""$minimiza" >> $arq_config
			echo "NumberProcessor = ""$nt" >> $arq_config
			echo "NumberObjective = ""$objetivos" >> $arq_config
			echo "NumberGeration = ""$geracoes" >> $arq_config
			echo "SizePopulation = ""$individuos" >> $arq_config
			echo "number_archive = ""$archive" >> $arq_config
			echo "NumberIndividualReproduce = ""$novosind" >> $arq_config
			echo "CrossoverRate = ""$crossover" >> $arq_config
			echo "blx_alfa = ""$blx_alpha" >> $arq_config
			echo "BLX_cros_Rate = ""$p_blx" >> $arq_config
			echo "1_point_cros_Rate = ""$p_1pt" >> $arq_config
			echo "2_point_cros_Rate = ""$p_2pt" >> $arq_config
			echo "MutationRate = ""$mutacao" >> $arq_config
			echo "Individual_Mutation_Rate = ""$mutacao_ind" >> $arq_config
			echo "max_mutation_range = ""$max_mut" >> $arq_config
	#		echo "rotamer_library = none" >> $arq_config
			echo "rotamer_library = cad_tuffery" >> $arq_config

			# Monta a string que contém os objetivos para ser colocada no arquivo de parâmetros
			i=1
			fitness_energy=""	# string começa vazia
			while [ $i -le $objetivos ]; do

				if [ $i -ne 1 ]; then
					fitness_energy="$fitness_energy"", "
				fi

				fitness_energy="$fitness_energy""${obj[$i]}"
				let i=$i+1
			done

			echo "Fitness_Energy = ""$fitness_energy" >> $arq_config

			echo "ArqFimMulti = ""$local""$titulo""_""$rep""/arqFim_""$proteina"".txt" >> $arq_config
			echo "NativeProtein = ""$local""$titulo""_""$rep""/""$proteina"".pdb" >> $arq_config
			echo "SequenceAminoAcidsPathFileName = ""$local""$titulo""_""$rep""/""$proteina"".fasta.txt" >> $arq_config
			echo "StoreFitnessResultsPathFileName = ""$local""$titulo""_""$rep""/result_""$proteina"".txt" >> $arq_config
			echo "PDBBestIndividualPathFileName = prot_""$proteina"".pdb" >> $arq_config >> $arq_config
			echo "FinalPopulationPathFileName = ""$local""$titulo""_""$rep""/""saida_""$proteina"".txt" >> $arq_config
			echo "NameExecutation = ""$titulo" >> $arq_config
			echo "Local_Execute = ""$local""$titulo""_""$rep""/" >> $arq_config
			echo "Database = ""$path_protpred""Database/" >> $arq_config
			echo "ComputeEnergyProgram = ""$path_protpred""scripts/compute_energy/run_gromacs_compute_energy.sh" >> $arq_config
			echo "MinimizationProgram = ""$path_protpred""scripts/compute_energy/run_gromacs_energy_minimization.sh" >> $arq_config
			echo "CleanGromacsSimulation = ""$path_protpred""scripts/compute_energy/clean_simulation.sh" >> $arq_config
			echo "GetEnergyProgram = ""$path_protpred""scripts/compute_energy/run_g_energy.sh" >> $arq_config
			echo "TopologyFile = topol_ProtPred.top" >> $arq_config
			echo "IniPopFileName = pop_""$proteina"".txt" >> $arq_config
			echo "z_matrix_fileName = z_matrix_1VII.z" >> $arq_config
			echo "Path_Gromacs_Programs = ""$path_gromacs" >> $arq_config
			echo "Computed_Energies_Gromacs_File = file_energy_computed.ener.edr" >> $arq_config
			echo "Energy_File_xvg = energy.xvg" >> $arq_config
			echo "Program_Read_Energy = ""$path_protpred""scripts/compute_energy/run_read_energy.sh" >> $arq_config
			echo "Program_Run_g_sas = ""$path_protpred""scripts/sas/run_g_sas.sh" >> $arq_config
			echo "GetAreasFrom_g_sas = ""$path_protpred""scripts/sas/run_read_areas.sh" >> $arq_config
			echo "Computed_Areas_g_sas_File = file_g_sas_areas.xvg" >> $arq_config
			echo "Computed_Energy_Value_File = energy_computed.txt" >> $arq_config
			echo "Program_Run_stride = ""$path_protpred""scripts/stride/run_stride.sh" >> $arq_config

			echo "objective_analisys =""$objective_analisys""" >> $arq_config
			echo "objective_analisys_dimo_source = ""$objective_analisys_dimo_source""" >> $arq_config	
			echo "Program_Run_GreedyTreeGenerator2PG = ""$path_protpred""/scripts/dimo/call_GreedyTreeGenerator2PG.sh" >> $arq_config



			# Monta string que tem os pesos dos objs para colocar no arquivo de parâmetros
			i=1
			weights=""
			while [ $i -le $objetivos ]; do
				if [ $i -ne 1 ]; then
					weights="$weights"", "
				fi

				weights="$weights""1.0"
				let i=$i+1
			done

			echo "Weights_Fitness = ""$weights" >> $arq_config
			echo "Program_Run_RMSD = ""$path_protpred""scripts/compute_rmsd/run_g_rms.sh" >> $arq_config
			echo "Program_Run_g_gyrate = ""$path_protpred""scripts/gyrate/run_g_gyrate.sh" >> $arq_config
			echo "GetRadiusFrom_g_gyrate = ""$path_protpred""scripts/gyrate/run_read_gyrate.sh" >> $arq_config
			echo "Computed_Radius_g_gyrate_File = file_g_gyrate_radius.xvg" >> $arq_config
			echo "Program_Run_g_hbond = ""$path_protpred""scripts/h_bond/run_g_hbond.sh" >> $arq_config
			echo "GetValueFrom_g_hbond = ""$path_protpred""scripts/h_bond/run_read_hbond.sh" >> $arq_config
			echo "Computed_g_hbond_File = file_g_hbond.xvg" >> $arq_config

		fi




		# Constroi o arquivo que contem os objetivos e o valor de maxmin que será necessário para alguns scripts de análise
		if [ -e "obj_maxmin" ]; then rm "obj_maxmin" ; fi	
		i=1
		while [ $i -le $objetivos ]; do	# para cada objetivo...

			if [ "${obj[$i]}" == "Potential" ]; then maxmin="min" ; fi
			if [ "${obj[$i]}" == "H_Bond_Main" ]; then maxmin="max" ; fi
			if [ "${obj[$i]}" == "H_Bond" ]; then maxmin="max" ; fi
			if [ "${obj[$i]}" == "Hydrophobic" ]; then maxmin="min" ; fi
			if [ "${obj[$i]}" == "Hydrophilic" ]; then maxmin="max" ; fi
			if [ "${obj[$i]}" == "Total_Area" ]; then maxmin="min" ; fi
			if [ "${obj[$i]}" == "Gyrate" ]; then maxmin="min" ; fi
			if [ "${obj[$i]}" == "Stride_helix" ]; then maxmin="max" ; fi
			if [ "${obj[$i]}" == "Stride_beta" ]; then maxmin="max" ; fi
			if [ "${obj[$i]}" == "Stride_total" ]; then maxmin="max" ; fi
			if [ "${obj[$i]}" == "GBSA_Solvatation" ]; then maxmin="min" ; fi



			echo "${obj[$i]} $maxmin" >> "obj_maxmin"

			let i=$i+1
		done





	# criação da população inicial
		if [ $POP_INI -eq 1 ]; then
			"$path_protpred""src/"./protpred-Gromacs_pop_initial $arq_config
		fi



	# algoritmo evolutivo
		if [ $EA -eq 1 ]; then

			if [ $objetivos -eq 1 ]; then
				"$path_protpred""src/"./protpred-Gromacs-Mono $arq_config
			fi


			if [ $objetivos -ge 2 ]; then
				if [ $algoritmo = "NSGA-II" ]; then
					"$path_protpred""src/"./protpred-Gromacs-NSGA2 $arq_config
				fi
			fi

		fi





	# ANÁLISES



	# análise mono
		if [ $MONO_ANALYSIS -eq 1 ]; then
			"$path_protpred""scripts/analysis/"./analysis_mono.sh "$path_protpred" "$arq_config" "$geracoes" "$individuos" "obj_maxmin"
		fi


	# análise Multi-Objetivo
		if [ $FRONT_ANALYSIS -eq 1 ]; then
			echo
			echo "Iniciando analise Multi-Objetivo - ""$titulo""_""$rep"
			"$path_protpred""scripts/analysis/"./fronts.sh "$path_protpred" "$path_gromacs" "$path_maxcluster" $individuos $geracoes $objetivos "$local""$titulo""_""$rep"/"$proteina"".pdb" "$arq_config" "obj_maxmin"
			echo "OK"
			echo
		fi

	# DM dos individuos nao-dominados da ultima geracao
		if [ $MOLECULAR_DYNAMICS -eq 1 ]; then
			echo
			echo "Iniciando simulacoes de dinamica molecular dos individuos nao-dominados - ""$titulo""_""$rep"
			"$path_protpred""scripts/analysis/"./md.sh "$path_protpred" "$path_gromacs" $individuos $geracoes "$arq_config" $nt "$local""$titulo""_""$rep"/"$proteina"".pdb"
			echo "OK"
			echo
		fi





	# sai da pasta da repeticao
		cd ../

		if [ $TAR -eq 1 ]; then
			if [ -e "$titulo""_""$rep"".tar.gz" ]; then rm "$titulo""_""$rep"".tar.gz" ; fi
		
			echo
			echo -n "Criando arquivo ""$titulo""_""$rep"".tar.gz""... "
			tar -czf "$titulo""_""$rep"".tar.gz" "$titulo""_""$rep""/"
			echo "OK"
		fi

		let rep=$rep+1

	done
# fim das repeticoes

	let pred_atual=$pred_atual+1

done
# fim das predicoes do arquivo de parametros
