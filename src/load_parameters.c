#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "load_parameters.h"
#include "LoadConfig.h"
#include "messages.h"
#include "string_owner.h"
#include "objective.h"


static void initialize_parameters_extended_chain(
		input_parameters_extended_chain_t *param){
	param->local_execute =  Malloc(char, MAX_PATH);
	param->seq_protein_file_name = Malloc(char, MAX_FILE_NAME);
	param->res_reference = -1;
	param->file_final_pdb = Malloc(char, MAX_FILE_NAME );
	param->path_database =  Malloc(char, MAX_PATH);
}

static void initialize_parameters(input_parameters_t *param){
	param->seq_protein_file_name = Malloc(char, MAX_PATH_FILE_NAME );
	param->file_final_pdb = Malloc(char, MAX_FILE_NAME );
	param->path_database =  Malloc(char, MAX_PATH );
	param->path_program_compute_energy = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_local_execute = Malloc(char, MAX_PATH );
	param->top_file = Malloc(char, MAX_FILE_NAME );
	param->initial_pop_file_name = Malloc(char, MAX_FILE_NAME );
	param->z_matrix_file = Malloc(char, MAX_FILE_NAME );
	param->path_program_minimization_program = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_gromacs_programs = Malloc(char, MAX_PATH );
	param->computed_energies_gromacs_file = Malloc(char, MAX_FILE_NAME );
	param->energy_file_xvg = Malloc(char, MAX_FILE_NAME );
	param->path_program_read_energy = Malloc(char, MAX_PATH_FILE_NAME );
	param->computed_energy_value_file = Malloc(char, MAX_FILE_NAME );
	param->computed_areas_g_sas_file = Malloc(char, MAX_FILE_NAME );
	param->computed_radius_g_gyrate_file = Malloc(char, MAX_FILE_NAME );
	param->computed_g_hbond_file = Malloc(char, MAX_FILE_NAME );
	param->path_program_g_energy = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_clean_simulation = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_file_native_protein = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_HIS_protonation = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_rmsd = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_read_xvg_rmsd = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_g_sas = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_read_g_sas = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_g_gyrate = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_read_g_gyrate = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_g_hbond = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_read_g_hbond = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_program_stride = Malloc(char, MAX_PATH_FILE_NAME );

	param->gromacs_energy_min = ener_min_none;
	param->gromacs_energy_min_gen_oper = ener_min_none;
	param->processor_number = 1;
	param->blx_alfa = 0.0;
	param->max_mutation_range = 0.010;
	param->individual_mutation_rate = 0.10;
	param->blx_cros_rate = -1;
	param->point_1_cros_rate = -1;
	param->point_2_cros_rate = -1;
	param->number_archive = -1;
	param->rotamer_library = rotamer_library_none;
	param->crossovers = NULL;
	param->objective_analysis = objective_analysis_none;
	param->path_dimo_sources = NULL;
	param->path_call_GreedyTreeGenerator2PG = NULL;
	param->number_individual_select_reproduce = 0;
}

static void set_number_individual_select_reproduce(input_parameters_t *param,
		const char *param_value){
	/*This function set the value for number_individual_select_reproduce
	 * It is based on number of immigrant individuals which is given by
	 * percentage_value.
	 * number_individual_select_reproduce parameter is the number of individuals
	 * which are selected to reproduce. The others will be added to population
	 * randomly. It will be the immigrant individuals.
	 */
	float perc = 0;
	int number_ind = 0;

	if (is_equal(param_value, "") == bfalse){
		perc = atof(param_value);
		number_ind = round((param->size_population*perc));
	}
	if (number_ind > param->size_population){
		char msg[1000];
		sprintf(msg, "In  set_number_individual_select_reproduce function, the number_ind is greater than size of population. Check the PercentageImmigrantIndivuals parameter must be between 0 and 1 \n",param_value);
		fatal_error(msg);
	}
	param->number_individual_select_reproduce = param->size_population - number_ind;
}

static void set_parameter_processor_number(input_parameters_t *param,
		const char *param_value){
	int aux;
	aux = atoi(param_value);
	if (aux > 1){
		param->processor_number = aux;
	}
}

static void set_parameter_gromacs_minimization(input_parameters_t *param,
		char *param_value){
	if (is_equal(param_value,"none") ){
		param->gromacs_energy_min = ener_min_none;
	}else if (is_equal(param_value,"ener_implicit") ){
		param->gromacs_energy_min = ener_min_implicit;
	}else if (is_equal(param_value,"ener_explicit") ){
		param->gromacs_energy_min = ener_min_explicit;
	}else{
		char msg[300];
		sprintf(msg, "gromacs_energy_min parameter was typed %s. But, this value must be none, ener_implicit or ener_explicit \n",param_value);
		fatal_error(msg);
	}
}

static void set_parameter_gromacs_minimization_gen_oper(input_parameters_t *param,
		char *param_value){
	if (is_equal(param_value,"none") ){
		param->gromacs_energy_min_gen_oper = ener_min_none;
	}else if (is_equal(param_value,"ener_implicit") ){
		param->gromacs_energy_min_gen_oper = ener_min_implicit;
	}else if (is_equal(param_value,"ener_explicit") ){
		param->gromacs_energy_min_gen_oper = ener_min_explicit;
	}else{
		char msg[300];
		sprintf(msg, "gromacs_energy_min_gen_oper parameter was typed %s. But, this value must be none, ener_implicit or ener_explicit \n",param_value);
		fatal_error(msg);
	}
}

static void set_parameter_fitness_energies(input_parameters_t *param,
		char *energy_parameters){
   /* Obtain the energies which will be used to compute the fitness
    * Receives parameter structure and the line contains the energies
    */
	int i = 0;
	char *c_params;
	param->fitness_energies = Malloc(type_fitness_energies_t
			,param->number_fitness);
	c_params = strtok(energy_parameters," ,");
	while (c_params != NULL){
		param->fitness_energies[i] = str2type_fitness_energies(c_params);
		i++;
		c_params = strtok(NULL," ,");
	}
	if (i != param->number_fitness){
		char msg[300];
		sprintf(msg, "The number of energies are %i. The number of fitness is %i. These numbers must be equal. Check it!", i , param->number_fitness);
		fatal_error(msg);
	}
}


static void set_parameter_weights_average(input_parameters_t *param,
		char *weights_parameters){
   /* Obtain the weights for fitness which will be used to compute the fitness
    * Receives parameter structure and the line contains the weights
    */
	int i = 0;
	char *c_params;
	param->weights_fitness = Malloc(float, param->number_fitness);
	c_params = strtok(weights_parameters," ,");
	while (c_params != NULL){
		param->weights_fitness[i] = str2float(c_params);
		i++;
		c_params = strtok(NULL," ,");
	}
	if (i != param->number_fitness){
		char msg[300];
		sprintf(msg, "The number of weights are %i. The number of fitness is %i. These numbers must be equal. Check it!", i , param->number_fitness);
		fatal_error(msg);
	}
}

static void set_parameter_rotamer_library(input_parameters_t *param,
		char *param_value){
	if (is_equal(param_value,"none") ){
		param->rotamer_library = rotamer_library_none;
	}else if (is_equal(param_value,"cad_tuffery") ){
		param->rotamer_library = rotamer_library_cad_tuffery;
	}else{
		char msg[300];
		sprintf(msg, "rotamer_library parameter was typed %s. But, this value must be either none or cad_tuffery \n",param_value);
		fatal_error(msg);
	}
}

static void set_parameter_objective_analysis(input_parameters_t *param,
		char *param_value, char *param_source_dimo,
		char *param_call_GreedyTreeGenerator2PG){
	if (is_equal(param_value,"none") ){
		param->objective_analysis = objective_analysis_none;
	}else if (is_equal(param_value,"dimoGreedyTreeGenerator") ){
		param->objective_analysis = objective_analysis_dimo_GreedyTreeGenerator;
		param->path_dimo_sources = Malloc(char, MAX_PATH_FILE_NAME );
		strcpy(param->path_dimo_sources, param_source_dimo);
		param->path_call_GreedyTreeGenerator2PG = Malloc(char, MAX_PATH_FILE_NAME );
		strcpy(param->path_call_GreedyTreeGenerator2PG, param_call_GreedyTreeGenerator2PG);
	}else{
		char msg[300];
		sprintf(msg, "objective_analisys parameter was typed %s. But, this value must be either none or dimoGreedyTreeGenerator \n",param_value);
		fatal_error(msg);
	}
}


void set_parameter_number_crossover(input_parameters_t *param){
	/* Computes how many crossovers were chosen to run the algorithm.
	*/
	if (param->crossover_rate > 0){ //means that crossover will be used.
		int how_many = 0;
		if (param->blx_cros_rate > 0){
			how_many = how_many +1;
		}
		if (param->point_1_cros_rate > 0){
			how_many = how_many +1;
		}
		if (param->point_2_cros_rate > 0){
			how_many = how_many +1;
		}
		if (how_many == 0){
			fatal_error("The number of crossover is 0. \n");
		}
		param->number_crossover = how_many;
	}
}

void set_crossovers(input_parameters_t *param){
	/* It stores the crossovers that the user chose.
	 */
	if (param->crossover_rate > 0){// It means that the crossover will be used.
		param->crossovers = Malloc(type_crossoers_t, param->number_crossover);

		int index = -1;
		if (param->blx_cros_rate > 0){
			index = index + 1;
			param->crossovers[index] = crossoer_blx_alpha;
		}
		if (param->point_1_cros_rate > 0){
			index = index + 1;
			param->crossovers[index] = crossoer_point_1;
		}
		if (param->point_2_cros_rate > 0){
			index = index + 1;
			param->crossovers[index] = crossoer_point_2;
		}
		if (index > param->number_crossover){
			fatal_error("In set_crossovers function there is an error of index. Check it!! \n");
		}
	}
}


void deAllocateload_parameters(input_parameters_t *param){
  
    free(param->seq_protein_file_name );
    free(param->file_final_pdb );
	free(param->path_database );
	free(param->path_program_compute_energy);
	free(param->path_local_execute);
	free(param->top_file );
	free(param->initial_pop_file_name );
	free(param->z_matrix_file);
	free(param->path_program_minimization_program);
	free(param->path_gromacs_programs);
	free(param->computed_energies_gromacs_file);
	free(param->energy_file_xvg);
	free(param->path_program_read_energy);
	free(param->computed_energy_value_file);
	free(param->computed_areas_g_sas_file);
	free(param->computed_radius_g_gyrate_file);
	free(param->computed_g_hbond_file);
	free(param->path_program_g_energy);
	free(param->path_program_clean_simulation);
	free(param->fitness_energies);
	free(param->weights_fitness);
	free(param->path_program_rmsd);
	free(param->path_program_read_xvg_rmsd);
	free(param->path_file_native_protein);
	free(param->path_program_HIS_protonation);
	free(param->path_program_g_sas);
	free(param->path_program_read_g_sas);
	free(param->path_program_g_gyrate);
	free(param->path_program_read_g_gyrate);
	free(param->path_program_g_hbond);
	free(param->path_program_read_g_hbond);
	free(param->path_program_stride);
	if (param->crossovers != NULL){
		free(param->crossovers);
	}
	if (param->path_dimo_sources != NULL){
		free(param->path_dimo_sources);
	}

	if (param->path_call_GreedyTreeGenerator2PG != NULL){
		free(param->path_call_GreedyTreeGenerator2PG);
	}
	//free(param);
}

void deAllocateload_parameters_extended_chain(
		input_parameters_extended_chain_t *param){
    free(param->seq_protein_file_name);
    free(param->file_final_pdb);
    free(param->path_database);
}

void load_parameters_from_file(input_parameters_t *param,
		const char *conf_file_name){
	/*Loading the configuration from file*/

	initialize_parameters(param);

	LoadConfig conf(conf_file_name);
	param->number_generation = atoi(conf.getParameterChar("NumberGeration"));
	param->size_population = atoi(conf.getParameter("SizePopulation").c_str());
	param->number_fitness = atoi(conf.getParameter("NumberObjective").c_str());
	param->crossover_rate = atof(conf.getParameter("CrossoverRate").c_str());
	param->mutation_rate = atof(conf.getParameter("MutationRate").c_str());
	strcpy(param->seq_protein_file_name, conf.getParameterChar("SequenceAminoAcidsPathFileName"));
	strcpy(param->file_final_pdb, conf.getParameterChar("PDBBestIndividualPathFileName"));
	strcpy(param->path_database, conf.getParameterChar("Database"));
	strcpy(param->path_program_compute_energy, conf.getParameterChar("ComputeEnergyProgram"));
	strcpy(param->path_local_execute, conf.getParameterChar("Local_Execute"));
	strcpy(param->top_file, conf.getParameterChar("TopologyFile"));
	strcpy(param->initial_pop_file_name, conf.getParameterChar("IniPopFileName"));
	strcpy(param->z_matrix_file, conf.getParameterChar("z_matrix_fileName"));
	strcpy(param->path_program_minimization_program, conf.getParameterChar("MinimizationProgram"));
	strcpy(param->path_gromacs_programs, conf.getParameterChar("Path_Gromacs_Programs"));
	strcpy(param->computed_energies_gromacs_file, conf.getParameterChar("Computed_Energies_Gromacs_File"));
	strcpy(param->energy_file_xvg, conf.getParameterChar("Energy_File_xvg"));
	strcpy(param->path_program_read_energy, conf.getParameterChar("Program_Read_Energy"));
	strcpy(param->computed_energy_value_file, conf.getParameterChar("Computed_Energy_Value_File"));
	strcpy(param->computed_areas_g_sas_file, conf.getParameterChar("Computed_Areas_g_sas_File"));
	strcpy(param->computed_radius_g_gyrate_file, conf.getParameterChar("Computed_Radius_g_gyrate_File"));
	strcpy(param->computed_g_hbond_file,conf.getParameterChar("Computed_g_hbond_File"));
	strcpy(param->path_program_g_energy,conf.getParameterChar("GetEnergyProgram"));
	strcpy(param->path_program_clean_simulation,conf.getParameterChar("CleanGromacsSimulation"));
	strcpy(param->path_file_native_protein,conf.getParameterChar("NativeProtein"));
	strcpy(param->path_program_HIS_protonation,conf.getParameterChar("HISProtonationGromacs"));
	strcpy(param->path_program_rmsd,conf.getParameterChar("Program_Run_RMSD"));
	strcpy(param->path_program_read_xvg_rmsd,conf.getParameterChar("Program_Read_RMSD"));
	strcpy(param->path_program_g_sas,conf.getParameterChar("Program_Run_g_sas"));
	strcpy(param->path_program_read_g_sas,conf.getParameterChar("GetAreasFrom_g_sas"));
	strcpy(param->path_program_g_gyrate,conf.getParameterChar("Program_Run_g_gyrate"));
	strcpy(param->path_program_read_g_gyrate,conf.getParameterChar("GetRadiusFrom_g_gyrate"));
	strcpy(param->path_program_g_hbond,conf.getParameterChar("Program_Run_g_hbond"));
	strcpy(param->path_program_read_g_hbond,conf.getParameterChar("GetValueFrom_g_hbond"));
	strcpy(param->path_program_stride,conf.getParameterChar("Program_Run_stride"));
	set_parameter_fitness_energies(param,conf.getParameterChar("Fitness_Energy"));
	set_parameter_weights_average(param,conf.getParameterChar("Weights_Fitness"));
	set_parameter_gromacs_minimization(param,conf.getParameterChar("gromacs_energy_min"));
	set_parameter_gromacs_minimization_gen_oper(param,conf.getParameterChar("gromacs_energy_min_gen_oper"));
	set_parameter_processor_number(param,conf.getParameter("NumberProcessor").c_str());
	param->blx_alfa = atof(conf.getParameter("blx_alfa").c_str());
	param->max_mutation_range = atof(conf.getParameter("max_mutation_range").c_str());
	param->individual_mutation_rate = atof(conf.getParameter("Individual_Mutation_Rate").c_str());
	param->blx_cros_rate = atof(conf.getParameter("BLX_cros_Rate").c_str());
	param->point_1_cros_rate = atof(conf.getParameter("1_point_cros_Rate").c_str());
	param->point_2_cros_rate = atof(conf.getParameter("2_point_cros_Rate").c_str());
	param->number_archive =  atoi(conf.getParameter("number_archive").c_str());

	set_parameter_rotamer_library(param,conf.getParameterChar("rotamer_library"));

	set_parameter_number_crossover(param);
	set_crossovers(param);

	set_parameter_objective_analysis(param,
			conf.getParameterChar("objective_analisys"),
			conf.getParameterChar("objective_analisys_dimo_source"),
			conf.getParameterChar("Program_Run_GreedyTreeGenerator2PG") );

	set_number_individual_select_reproduce(param,conf.getParameterChar("PercentageImmigrantIndivuals"));
	
}

void set_parameters_extended_chain(input_parameters_extended_chain_t *param,
		char *argv[]){
	initialize_parameters_extended_chain(param);

	/* Setting values from terminal */
	strcpy(param->local_execute, argv[1]);
	strcpy(param->seq_protein_file_name, argv[2]);
	strcpy(param->file_final_pdb, argv[3]);
	param->res_reference = str2int(argv[4]);
	strcpy(param->path_database, argv[5]);
	param->steps = str2int(argv[6]);
}

