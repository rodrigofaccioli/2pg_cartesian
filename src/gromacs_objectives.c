#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "gromacs_objectives.h"
#include "defines.h"
#include "messages.h"
#include "futil.h"
#include "consts.h"
#include "string_owner.h"
#include "osutil.h"
#include "parameters_type.h"
#include "protein_type.h"
#include "pdbio.h"


#define TAM_LINE_ENER 50
#define MAX_ENERGY 999999999.99
#define MIN_ENERGY -999999999.99

static char gyrate_file[MAX_FILE_NAME];

//It is based on GROMACS version 4.5.3
static option_fitness_gromacs_t option_g_energy_program [] = {
      		                                                  {gmx_potential_ener, "12","Potential"},
      		                                                  {gmx_edw_ener,"8","LJ-14"},
      		                                                  {gmx_elel_ener,"9","Coulomb-14"},
      		                                                  {gmx_hydrophobic,"-1","Hydrophobic"},
      		                                                  {gmx_hydrophilic,"-1","Hydrophilic"},
      		                                                  {gmx_total_area,"-1","Total_Area"},
      		                                                  {gmx_gyrate,"-1","Gyrate"},
      		                                                  {gmx_hbond,"-1","H_Bond"},
      		                                                  {gmx_hbond_main,"-1","H_Bond_Main"},
      		                                                  {gmx_GBSA_Solvatation,"-1","GBSA_Sol"},
      		                                                  {gmx_GB_Polarization,"6","GB-Polarization"},
      		                                                  {gmx_Nonpolar_Sol,"7","Nonpolar-Sol."},
      		                                                  {gmx_stride_total,"-1","Stride_total"},
      		                                                  {gmx_stride_helix,"-1","Stride_helix"},
      		                                                  {gmx_stride_beta,"-1","Stride_beta"}
                                                             };


option_g_energy_t get_option_g_energy_t_from_type_fitness_energy(const type_fitness_energies_t *fit_ener){
	/*Receives type_fitness_energies_t and returns its option_g_energy_t
	 * For example: fit_ener_potential is gmx_potential_ener
	 */
	if (*fit_ener == fit_ener_potential){
		return gmx_potential_ener;
	}else if (*fit_ener == fit_ener_edw){
		return gmx_edw_ener;
	}else if (*fit_ener == fit_ener_ele){
		return gmx_elel_ener;
	}else if (*fit_ener == fit_hydrophobic){
		return gmx_hydrophobic;
	}else if (*fit_ener == fit_hydrophilic){
		return gmx_hydrophilic;
	}else if (*fit_ener == fit_total_area){
		return gmx_total_area;
	}else if (*fit_ener == fit_gyrate){
		return gmx_gyrate;
	}else if (*fit_ener == fit_hbond){
		return gmx_hbond;
	}else if (*fit_ener == fit_hbond_main){
		return gmx_hbond_main;
	}else if (*fit_ener == fit_GBSA_Solvatation){
		return gmx_GBSA_Solvatation;
	}else if (*fit_ener == fit_stride_total){
		return gmx_stride_total;
	}else if (*fit_ener == fit_stride_helix){
		return gmx_stride_helix;
	}else if (*fit_ener == fit_stride_beta){
		return gmx_stride_beta;
	}else{
		fatal_error("Option did not find at get_option_g_energy_t_from_type_fitness_energy function. Please check it! ");
	}
}

/* Copies to command the program which cleans the simulation.
*/
static void build_clean_command(char * command,
		const char *path_clean_simulation_program){
	strcpy(command,path_clean_simulation_program);
	strcat(command, " > /dev/null 2> /dev/null");
}

/** Calls script which remove files used in simulation.
* Simulation means an execution to calculate the objectives by GROMACS.
*/
static int clean_gromacs_simulation(const char *path_clean_simulation_program){
	int ret;
	char *command;

	command = Malloc(char,MAX_COMMAND);
	build_clean_command(command, path_clean_simulation_program);

	ret = system(command);
    free(command);

	return ret;
}

/* Sets the value of objective in solution 
* sol represents the solution (individual). The value of objective will be 
*     set in obj_values field.
* path_file_ener is the name of file that contain the value of objective
* fit represents the index of objective
*/
static void set_objective_from_gromacs_file_in_solution(solution_t *sol,char *path_file_ener,
		const int *fit){
	char *line;
	line = Malloc(char, TAM_LINE_ENER);
	strcpy(line, " ");
	double aux = -1;
	FILE *file_ener = open_file(path_file_ener,fREAD);
	if (*fit > sol->num_obj){
		fatal_error("Error when try to execute set_objective_from_gromacs_file_in_solution function \n");
	}
	fgets(line, TAM_LINE_ENER, file_ener);
	//sscanf(line,"%f",&aux);
	trim(line);
	aux = str2double(line);
	if (aux ==  -1.000000){ //This value means that Gromacs considered the energy as an infinitive value
		aux = MAX_ENERGY;
	}
	sol->obj_values[*fit] = aux;
	fclose(file_ener);
	free(line);
}

static void build_gromacs_command_for_gyrate(char *command,
		const char *local_execute, const char *computed_g_gyrate_program,
		const char *path_gromacs_programs, const char *pdbfile,
		const char *computed_g_gyrate_value_file){
	char space [] = " ";

	strcpy(command,computed_g_gyrate_program);
	strcat(command,space);
	strcat(command,path_gromacs_programs);
	strcat(command,space);
	strcat(command,local_execute);
	strcat(command,space);
	strcat(command,pdbfile);
	strcat(command,space);
	strcat(command,computed_g_gyrate_value_file);
}

/** Creates file name that received from xvg file the value of radiius of protein
 */
static void build_file_gyrate(){
	strcpy(gyrate_file,option_g_energy_program[gmx_gyrate].option_name);
	strcat(gyrate_file,".radii");
}

/** Builds command line for calling script which read a xvg file
* and saving radius of protein from g_gyrate program
*/
static void build_command_read_xvg_gyrate_file(char *command,
				const char *local_execute, const char *path_program_read_g_gyrate,
				const char *computed_g_gyrate_value_file){
	char space [] = " ";

	strcpy(command,path_program_read_g_gyrate);
	strcat(command,space);
	strcat(command,local_execute);
	strcat(command,space);
	strcat(command,computed_g_gyrate_value_file); //xvg file
	strcat(command,space);
	strcat(command,gyrate_file);
}

		
/** runs g_gyrate program and evaluates the radius of gyration of the protein
 * based on alpha-carbons only
 */
int compute_gyrate_mono(solution_t *sol, const char *local_execute,
		const char *path_gromacs_programs, const char *pdbfile,
		const option_g_energy *opt_fitness,
		const char *computed_g_gyrate_program,
		const char *path_program_read_g_gyrate,
		const char *computed_g_gyrate_value_file,
		const char *path_program_clean_simulation){
	int ret;
	int fit = 0;
	char *command;
	char *gyrate_path_file_name;

	command = Malloc(char,MAX_COMMAND);
	build_gromacs_command_for_gyrate(command, local_execute, computed_g_gyrate_program,
			path_gromacs_programs, pdbfile, computed_g_gyrate_value_file);
	ret = system(command); //calls to run g_gyrate
	if (check_exists_file(computed_g_gyrate_value_file) == btrue){
		build_file_gyrate();
		build_command_read_xvg_gyrate_file(command,local_execute,
				path_program_read_g_gyrate, computed_g_gyrate_value_file);
		ret = system(command); //calls to read g_gyrate xvg
		gyrate_path_file_name = path_join_file(local_execute,gyrate_file);
		//Set the value of radius option in solution
	    set_objective_from_gromacs_file_in_solution(sol,gyrate_path_file_name, &fit);
	    free(gyrate_path_file_name);

	}else{
		sol->obj_values[fit] = MAX_ENERGY; //MAX energy value
	}

	clean_gromacs_simulation(path_program_clean_simulation);

	free(command);
}

static void build_pdb_file_name(char *pdb_file_name, const char *aux_name,
		const char *__restrict prefix){
	strcpy(pdb_file_name, prefix);
	strcat(pdb_file_name,aux_name);
	strcat(pdb_file_name,".pdb");
}

/** Calculates the objectives by GROMACS
*/
void get_gromacs_objectives(solution_t *solutions, const input_parameters_t *in_para){
	const protein_t *population_aux;
	char pdbfile_aux[30];
	char aux [20];
	char aux_ind[5];
	char msg [50];	
	option_g_energy *opt_objective;

	//getting the gromacs objectives from input parameters
	opt_objective = Malloc(option_g_energy,in_para->number_fitness);
	for (int ob = 0; ob < in_para->number_fitness; ob++){
		opt_objective[ob] = get_option_g_energy_t_from_type_fitness_energy(&in_para->fitness_energies[ob]);
	}
	// calculating the objectivies of population
	for (int ind = 0; ind < in_para->size_population; ind++){
		//getting the protein from solution
		population_aux = (protein_t*) solutions[ind].representation;
		// saving pdb file of protein
    	int2str(aux_ind,&ind);
    	build_pdb_file_name(pdbfile_aux, aux_ind, PREFIX_PDB_FILE_NAME_EA);
   	    save_pdb_file(in_para->path_local_execute, pdbfile_aux, 
   	    	&population_aux->p_topol->numatom, population_aux->p_atoms, NULL);
   	    // obtaing the values of objectivies
   	    for (int ob = 0; ob < in_para->number_fitness; ob++){
	   	    if (opt_objective[ob] == gmx_gyrate) {
   		    	compute_gyrate_mono(&solutions[ind],in_para->path_local_execute,
   	    				in_para->path_gromacs_programs,pdbfile_aux,opt_objective,
   	    				in_para->path_program_g_gyrate,
   	    				in_para->path_program_read_g_gyrate,
   	    				in_para->computed_radius_g_gyrate_file,
   	    				in_para->path_program_clean_simulation);
   	    	}
   	    }
	}

	free(opt_objective);

}