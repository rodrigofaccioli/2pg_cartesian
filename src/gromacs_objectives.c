#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <errno.h>
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
#define MAX_ENERGY 999999999999999999999999999999.9999
#define MIN_ENERGY -99999999999999999999999999999.9999
#define MAX_VALUE 10

static char *program = NULL; /* program (executable) file name */
static char *filenm1 = NULL; /* data file names */
static char *filenm2 = NULL;
static char *filenm3 = NULL;
static char *filenm4 = NULL;
static char *filenm5 = NULL;


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


static inline int run_program(const char *file, char *const argv[])
{
	pid_t pid;
	int status;
	int out, err; /* file descriptors for stdout and stderr */

	pid = fork();

	if (pid == -1) {
		perror("Fork failed to create process");
		return 0;
	} else if (pid == 0) {
		/* child process */
		out = open("/dev/null", O_RDONLY);
		err = open("/dev/null", O_RDONLY);
		dup2(out, 1);
		dup2(err, 2);

		execv(file, argv);
		perror("Execv failed to run program (run_program)");
		_exit(EXIT_FAILURE);
	} else {
		/* parent process */
		waitpid(pid, &status, 0);
	}

	return 1;
}

/* Start a program as in:
 * echo "pipe_msg" | program [args...] */
static inline int run_program_after_pipe(const char *pipe_msg, const char *file,
							char *const argv[])
{
	pid_t pid;
	int status;
	int fd[2], out, err;

	if (pipe(fd) == -1) {
		perror("Failed to create pipe");
		return 0;
	}

	pid = fork();

	if (pid == -1) {
		perror("Fork failed to create process");
		return 0;
	} else if (pid == 0) {
		/* child process */
		/* supress output */
		out = open("/dev/null", O_WRONLY);
		err = open("/dev/null", O_WRONLY);
		dup2(out, 1);
		dup2(err, 2);

		/* read from the pipe */
		close(fd[1]);
		dup2(fd[0], 0);
		execv(file, argv);
		perror("Execv failed to run program");
		_exit(EXIT_FAILURE);
	} else {
		/* parent process */
		close(fd[0]);
		write(fd[1], pipe_msg, strlen(pipe_msg) + 1);
		close(fd[1]);
		waitpid(pid, &status, 0);
	}

	return 1;
}

/*
 * Run nprogs programs with pipes interconecting them and (optionally, if
 * output_file is not NULL) write the output to a file, as in:
 * 
 * prog0 [args0] | prog1 [args1] | ... | progn-1 [argsn-1] > output_file
 *
 * where:
 * argv_list[i] = argsi (argsi is in the format expected by execv(3))
 */
static inline int run_programs_with_pipe(int nprogs, char ***const argv_list,
							const char *output_file)
{
	int oldpipe[2], newpipe[2];
	int status;
	int i;
	int ret_value;
	char **argv;

	if (pipe(oldpipe) == -1) {
		perror("Failed to create pipe");
		return 0;
	}

	ret_value = 1;
	/* create first child */
	switch(fork()) {
		case -1:
			perror("Fork failed to create process");
			return 0;

		case 0:
			/* 1st child process writes to oldpipe, reads from none */
			close(oldpipe[0]);
			dup2(oldpipe[1], 1);

			argv = argv_list[0];
			execv(argv[0], argv);
			perror("Execv failed to run program");
			_exit(EXIT_FAILURE);
	}

	/* create intermediate children */
	/* 1st child process already created, last one to be created later */
	for (i = 1; i < nprogs-1; i++) {
		if (pipe(newpipe) == -1) {
			perror("Failed to create pipe");
			ret_value = 0;
			goto end;
		}

		switch (fork()) {
			case -1:
				perror("Fork failed to create process");
				ret_value = 0;
				goto end;

			case 0:
				/* child reads from oldpipe, writes to newpipe */
				close(oldpipe[1]);
				close(newpipe[0]);
				dup2(oldpipe[0], 0);
				dup2(newpipe[1], 1);

				argv = argv_list[i];
				execv(argv[0], argv);
				perror("Execv failed to run program");
				_exit(EXIT_FAILURE);

			default:
				/* parent process */
				close(oldpipe[0]);
				close(oldpipe[1]);
				oldpipe[0] = newpipe[0];
				oldpipe[1] = newpipe[1];
		}
	}

	/* create last child */
	switch(fork()) {
		case -1:
			perror("Fork failed to create process");
			ret_value = 0;
			goto end;

		case 0:
			/* last child process reads from oldpipe, maybe writes to file */
			close(oldpipe[1]);
			dup2(oldpipe[0], 0);

			if (output_file != NULL) {
				int outfd;

				outfd = open(output_file, O_WRONLY | O_CREAT | O_TRUNC,
						S_IRUSR | S_IWUSR | S_IRGRP);
				if (outfd == -1) {
					perror("Error opening output file");
					_exit(EXIT_FAILURE);
				}
				dup2(outfd, 1);
			}

			argv = argv_list[nprogs-1];
			execv(argv[0], argv);
			perror("Execv failed to run program");
			_exit(EXIT_FAILURE);
	}

	close(oldpipe[0]);
	close(oldpipe[1]);

	end:
	for (; i >= 0; i--) /* wait for all the created children processes */
		wait(&status);

	return ret_value;
}


/** Initialize GROMACS execution
* This function must be called exactly once before any other one in this file
 * to obtain memory for internal variables used to interact with Gromacs.
 *
 * Call finish_gromacs_execution when done using Gromacs
 */
void init_gromacs_execution (){
	program = Malloc(char, MAX_COMMAND);
	filenm1 = Malloc(char, MAX_COMMAND);
	filenm2 = Malloc(char, MAX_COMMAND);
	filenm3 = Malloc(char, MAX_COMMAND);
	filenm4 = Malloc(char, MAX_COMMAND);
	filenm5 = Malloc(char, MAX_COMMAND);
}

/** Finish GROMACS execution 
* This function is to be called when the Gromacs funcions below will no longer
* be used in order to release memory used by internal variables 
*/
void finish_gromacs_execution(){
	free(program);
	free(filenm1);
	free(filenm2);
	free(filenm3);
	free(filenm4);
	free(filenm5);

	program = NULL;
	filenm1 = NULL;
	filenm2 = NULL;
	filenm3 = NULL;
	filenm4 = NULL;
	filenm5 = NULL;
}


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
* value the calculated value that will be set in solution
* fit represents the index of objective
*/
static void set_objective_from_gromacs_in_solution(solution_t *sol,char *value,
		const int *fit){
	double aux = -1;

	aux = str2double(value);
	if (aux ==  -1.000000){ //This value means that Gromacs considered the energy as an infinitive value
		aux = MAX_ENERGY;
	}
	sol->obj_values[*fit] = aux;
	
	
}


/** runs g_gyrate program and evaluates the radius of gyration of the protein
 * based on alpha-carbons only
 */
static int compute_gyrate(solution_t *sol, const int *fit, const char *local_execute,
		const char *path_gromacs_programs, const char *pdbfile,
		const option_g_energy *opt_fitness, const char *computed_g_gyrate_value_file,
		const char *path_program_clean_simulation){
	int ret;
	
	char *last_line, *line_splited;
	char *value;
	char *opt_f;
	char *opt_s;
	char *opt_o;

	char *g_gyrate_args[8];
	opt_f = Malloc(char,2);
	opt_s = Malloc(char,2);
	opt_o = Malloc(char,2);
	
	value = Malloc(char,MAX_VALUE);

	//Setting variables to use
	strcpy(opt_f, "-f");
	strcpy(opt_s, "-s");
	strcpy(opt_o, "-o");

	/* g_gyrate - radius of gyration */
	strcpy(program, path_gromacs_programs);
	strcat(program, "g_gyrate");
	g_gyrate_args[0] = program;
	//pdb
	g_gyrate_args[1] = opt_f;
	strcpy(filenm1, local_execute);
	strcat(filenm1, pdbfile);
	g_gyrate_args[2] = filenm1;
	//tpr
	g_gyrate_args[3] = opt_s;
	strcpy(filenm2, local_execute);
	strcat(filenm2, "prot.tpr"); /* local_execute/prot.tpr */
	g_gyrate_args[4] = filenm2;
	//xvg
	g_gyrate_args[5] = opt_o;
	strcpy(filenm3, local_execute);
	strcat(filenm3, computed_g_gyrate_value_file);	
	g_gyrate_args[6] = filenm3;
	//final parameter
	g_gyrate_args[7] = NULL;

	if (!run_program_after_pipe("C-alpha", program, g_gyrate_args))
		fatal_error("Failed to run g_gyrate\n");

	if (check_exists_file(computed_g_gyrate_value_file) == btrue){
		//get the last line of xvg file
		last_line = get_last_line(computed_g_gyrate_value_file);
		// Split last line by space and obtaing the first value
 		line_splited = strtok (last_line," ");
 		// Obtaining the second value. It will be set in solution
 		line_splited = strtok(NULL, " ");
 		strcpy(value, line_splited);
 		//Looking the end of line_splited
	  	while (line_splited != NULL){	    	
    		line_splited = strtok(NULL, " ");
  		}
		//Set the value of radius option in solution
	    set_objective_from_gromacs_in_solution(sol,value, fit);
	    free(line_splited);
	    free(last_line);
	    free(value);
	}else{
		sol->obj_values[*fit] = MAX_ENERGY; //MAX energy value
	}

	clean_gromacs_simulation(path_program_clean_simulation);
	free(opt_f);
	free(opt_s);
	free(opt_o);
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
   		    	compute_gyrate(&solutions[ind], &ob, in_para->path_local_execute,
   	    				in_para->path_gromacs_programs, pdbfile_aux, opt_objective,
   	    				in_para->computed_radius_g_gyrate_file,
   	    				in_para->path_program_clean_simulation);
   	    	}
   	    }
	}

	free(opt_objective);
}