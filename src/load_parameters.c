#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "load_parameters.h"
#include "LoadConfig.h"
#include "messages.h"
#include "string_owner.h"
#include "objective.h"


static void initialize_parameters(input_parameters_t *param){
	param->seq_protein_file_name = Malloc(char, MAX_PATH_FILE_NAME );
	param->path_local_execute = Malloc(char, MAX_PATH );
	param->initial_pop_file_name = Malloc(char, MAX_FILE_NAME );
	param->path_gromacs_programs = Malloc(char, MAX_PATH );
	param->computed_energies_gromacs_file = Malloc(char, MAX_FILE_NAME );
	param->energy_file_xvg = Malloc(char, MAX_FILE_NAME );
	param->computed_energy_value_file = Malloc(char, MAX_FILE_NAME );
	param->computed_areas_g_sas_file = Malloc(char, MAX_FILE_NAME );
	param->computed_radius_g_gyrate_file = Malloc(char, MAX_FILE_NAME );
	param->computed_g_hbond_file = Malloc(char, MAX_FILE_NAME );
	param->individual_mutation_rate = 0.10;
	param->objective_analysis = objective_analysis_none;
	param->path_dimo_sources = NULL;
	param->path_call_GreedyTreeGenerator2PG = NULL;	

	param->crossovers = NULL;
	
    param->min_angle_mutation_phi = -180.0;
	param->max_angle_mutation_phi = 180.0;
    
    param->min_angle_mutation_psi = -180.0;
    param->max_angle_mutation_psi = 180.0;
    
    param->min_angle_mutation_omega = -180.0;
    param->max_angle_mutation_omega = 180.0;
    
    param->min_angle_mutation_side_chain = -180.0;
    param->max_angle_mutation_side_chain = 180.0;

    param->force_field = Malloc(char, MAX_FORCE_FIELD_NAME);    
    param->mdp_file = Malloc(char, MAX_FILE_NAME);
    
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
		param->fitness_energies[i] = str2type_objective(c_params);
		i++;
		c_params = strtok(NULL," ,");
	}
	if (i != param->number_fitness){
		char msg[300];
		sprintf(msg, "The number of energies are %i. The number of fitness is %i. These numbers must be equal. Check it!", i , param->number_fitness);
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
	int how_many = 1;
	if (how_many == 0){
		fatal_error("The number of crossover is 0. \n");
	}			
	param->number_crossover = how_many;
}

void set_crossovers(input_parameters_t *param){
	/* It stores the crossovers that the user chose.
	 */
	/*
	if (param->crossover_rate > 0){// It means that the crossover will be used.
		param->crossovers = Malloc(type_crossoers_t, param->number_crossover);

		int index = -1;
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
	*/
}


void deAllocateload_parameters(input_parameters_t *param){
    free(param->seq_protein_file_name );
	free(param->initial_pop_file_name );
	free(param->computed_energies_gromacs_file);
	free(param->energy_file_xvg);
	free(param->computed_energy_value_file);
	free(param->computed_areas_g_sas_file);
	free(param->computed_radius_g_gyrate_file);
	free(param->computed_g_hbond_file);
	free(param->fitness_energies);
	free(param->force_field);
	free(param->mdp_file);
	if (param->crossovers != NULL){
		free(param->crossovers);
	}
	if (param->path_dimo_sources != NULL){
		free(param->path_dimo_sources);
	}

	if (param->path_call_GreedyTreeGenerator2PG != NULL){
		free(param->path_call_GreedyTreeGenerator2PG);
	}
}

void load_parameters_from_file(input_parameters_t *param,
		const char *conf_file_name){
	/*Loading the configuration from file*/

	initialize_parameters(param);

	LoadConfig conf(conf_file_name);	
	param->number_generation = atoi(conf.getParameterChar("NumberGeration"));
	param->size_population = atoi(conf.getParameter("SizePopulation").c_str());
	param->number_fitness = atoi(conf.getParameter("NumberObjective").c_str());	
	strcpy(param->path_gromacs_programs, conf.getParameterChar("Path_Gromacs_Programs"));
	strcpy(param->seq_protein_file_name, conf.getParameterChar("SequenceAminoAcidsPathFileName"));	
	strcpy(param->path_local_execute, conf.getParameterChar("Local_Execute"));
	strcpy(param->initial_pop_file_name, conf.getParameterChar("IniPopFileName"));
	strcpy(param->computed_energies_gromacs_file, conf.getParameterChar("Computed_Energies_Gromacs_File"));
	strcpy(param->energy_file_xvg, conf.getParameterChar("Energy_File_xvg"));
	strcpy(param->computed_energy_value_file, conf.getParameterChar("Computed_Energy_Value_File"));
	strcpy(param->computed_areas_g_sas_file, conf.getParameterChar("Computed_Areas_g_sas_File"));
	strcpy(param->computed_radius_g_gyrate_file, conf.getParameterChar("Computed_Radius_g_gyrate_File"));
	strcpy(param->computed_g_hbond_file,conf.getParameterChar("Computed_g_hbond_File"));
	set_parameter_fitness_energies(param,conf.getParameterChar("Fitness_Energy"));
	param->individual_mutation_rate = atof(conf.getParameter("Individual_Mutation_Rate").c_str());
	param->point_1_cros_rate = atof(conf.getParameter("1_point_cros_Rate").c_str());
	param->point_2_cros_rate = atof(conf.getParameter("2_point_cros_Rate").c_str());

	set_parameter_number_crossover(param);
	set_crossovers(param);

	set_parameter_objective_analysis(param,
			conf.getParameterChar("objective_analisys"),
			conf.getParameterChar("objective_analisys_dimo_source"),
			conf.getParameterChar("Program_Run_GreedyTreeGenerator2PG") );


    param->min_angle_mutation_phi = atof(conf.getParameter("min_angle_mutation_phi").c_str());
	param->max_angle_mutation_phi = atof(conf.getParameter("max_angle_mutation_phi").c_str());

    param->min_angle_mutation_psi = atof(conf.getParameter("min_angle_mutation_psi").c_str());
    param->max_angle_mutation_psi = atof(conf.getParameter("max_angle_mutation_psi").c_str());

    param->min_angle_mutation_omega = atof(conf.getParameter("min_angle_mutation_omega").c_str());
    param->max_angle_mutation_omega = atof(conf.getParameter("max_angle_mutation_omega").c_str());

    param->min_angle_mutation_side_chain = atof(conf.getParameter("min_angle_mutation_side_chain").c_str());
    param->max_angle_mutation_side_chain = atof(conf.getParameter("max_angle_mutation_side_chain").c_str());


	strcpy(param->mdp_file,conf.getParameterChar("mdp_file_name"));
	strcpy(param->force_field,conf.getParameterChar("force_field"));

}