#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#ifndef WIN32
#include <sys/wait.h>
#include <unistd.h>
#else
#include <float.h>
#endif

#include <fcntl.h>
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

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif


#define TAM_LINE_ENER 50
#define MAX_ENERGY 9999999999999999999999999.9999				   
#define MIN_ENERGY -9999999999999999999999999.9999
#define MAX_VALUE 40
#define MAX_VALUE_G_SAS 3

static char *command = NULL;
static char *program = NULL; /* program (executable) file name */
static char *filenm1 = NULL; /* data file names */
static char *filenm2 = NULL;
static char *filenm3 = NULL;
static char *filenm4 = NULL;
static char *filenm5 = NULL;
//Options for GROMACS programs
static char *f_step0 = NULL;
static char *prot_gro = NULL;
static char *prot_top = NULL;
static char *prot_tpr = NULL;
static char *prot_trr = NULL;
static char *prot_log = NULL;
static char *confout_gro = NULL;
static char *posre_itp = NULL;
static char *mdout_mdp = NULL;
static char *file_energy_computed_ener_edr = NULL;
static char *prot_sys_trr = NULL;
static char *prot_sys_tpr = NULL;
static char *prot_sys_top = NULL;
static char *prot_sys_gro = NULL;
static char *energy_xvg = NULL;
static char *traj_xtc = NULL;
static char *opt_f = NULL;
static char *opt_s = NULL;
static char *opt_o = NULL;
static char *opt_c = NULL;
static char *opt_ff = NULL;
static char *opt_water = NULL;
static char *opt_none = NULL;
static char *opt_p = NULL;
static char *opt_ignh = NULL;
static char *xvg_1 = NULL;
static char *opt_rerun = NULL;
static char *opt_e = NULL;
static char *opt_g = NULL;
static char *opt_num = NULL;
/* Stores the values of g_sas program 
* g_sas_values[0] Hydrophobic
* g_sas_values[1] Hydrophilic
* g_sas_values[3] Total Area
*/
static double *g_sas_values = NULL;


//It is based on GROMACS version 4.6.5
static option_fitness_gromacs_t option_g_energy_program [] = {
      		                                                  {gmx_potential_ener, "11","Potential"},
      		                                                  {gmx_edw_ener,"7","LJ-14"},
      		                                                  {gmx_elel_ener,"8","Coulomb-14"},
      		                                                  {gmx_hydrophobic,"-1","Hydrophobic"},
      		                                                  {gmx_hydrophilic,"-1","Hydrophilic"},
      		                                                  {gmx_total_area,"-1","Total_Area"},
      		                                                  {gmx_gyrate,"-1","Gyrate"},
      		                                                  {gmx_hbond,"-1","H_Bond"},
      		                                                  {gmx_hbond_main,"-1","H_Bond_Main"},
      		                                                  {gmx_GBSA_Solvatation,"-1","GBSA_Sol"},
      		                                                  {gmx_GB_Polarization,"5","GB-Polarization"},
      		                                                  {gmx_Nonpolar_Sol,"6","Nonpolar-Sol."},
      		                                                  {gmx_stride_total,"-1","Stride_total"},
      		                                                  {gmx_stride_helix,"-1","Stride_helix"},
      		                                                  {gmx_stride_beta,"-1","Stride_beta"}
                                                             };


static inline int run_program(const char *file, char *const argv[]){
	int out, err; /* file descriptors for stdout and stderr */
	int i = 0;

	strcpy(command, argv[i]);
	strcat(command, " ");
	i++;
	do{	
		strcat(command, argv[i]);
		strcat(command, " ");
		i++;	
	}while (argv[i] != NULL);
	//Avoid output messages
	strcat(command, " > /dev/null 2> /dev/null ");

	system(command);

	return 1;
}

/* Start a program as in:
 * echo "pipe_msg" | program [args...] */
static inline int run_program_after_pipe(const char *pipe_msg, const char *file,
							char *const argv[]){
	int out, err; /* file descriptors for stdout and stderr */
	int i;

	strcpy(command, "echo ");
	strcat(command, pipe_msg);
	strcat(command, " | ");

	i = 0;
	do{	
		strcat(command, argv[i]);
		strcat(command, " ");
		i++;	
	}while (argv[i] != NULL);
	//Avoid output messages
	strcat(command, " > /dev/null 2> /dev/null ");

	system(command);
	
	return 1;
}

#ifndef WIN32
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
#endif

void initialize_g_sas_values(){
	//Hydrophobic
	g_sas_values[0] = MAX_ENERGY;
	//Hydrophilic
	g_sas_values[1] = MIN_ENERGY;
	//Total Area 
	g_sas_values[2] = MAX_ENERGY;
}

/** Initialize GROMACS execution
* This function must be called exactly once before any other one in this file
 * to obtain memory for internal variables used to interact with Gromacs.
 *
 * Call finish_gromacs_execution when done using Gromacs
 */
_2PG_CARTESIAN_EXPORT
void init_gromacs_execution (){

	command = Malloc(char, MAX_COMMAND);
	program = Malloc(char, MAX_COMMAND);
	filenm1 = Malloc(char, MAX_COMMAND);
	filenm2 = Malloc(char, MAX_COMMAND);
	filenm3 = Malloc(char, MAX_COMMAND);
	filenm4 = Malloc(char, MAX_COMMAND);
	filenm5 = Malloc(char, MAX_COMMAND);

	f_step0 = Malloc(char, MAX_FILE_NAME);
	prot_gro = Malloc(char, MAX_FILE_NAME);
	prot_top = Malloc(char, MAX_FILE_NAME);
	prot_tpr = Malloc(char, MAX_FILE_NAME);
	prot_trr = Malloc(char, MAX_FILE_NAME);
	prot_log = Malloc(char, MAX_FILE_NAME);
	confout_gro = Malloc(char, MAX_FILE_NAME);
	posre_itp = Malloc(char, MAX_FILE_NAME);
	mdout_mdp = Malloc(char, MAX_FILE_NAME);
	file_energy_computed_ener_edr = Malloc(char, MAX_FILE_NAME);
	prot_sys_trr = Malloc(char, MAX_FILE_NAME);
	prot_sys_tpr = Malloc(char, MAX_FILE_NAME);
	prot_sys_top = Malloc(char, MAX_FILE_NAME);
	prot_sys_gro = Malloc(char, MAX_FILE_NAME);
	energy_xvg = Malloc(char, MAX_FILE_NAME);
	traj_xtc = Malloc(char, MAX_FILE_NAME);
	xvg_1  = Malloc(char, MAX_FILE_NAME);
	opt_f = Malloc(char,3);
	opt_s = Malloc(char,3);
	opt_o = Malloc(char,3);
	opt_ff  = Malloc(char, 4);
	opt_water = Malloc(char, 7);
	opt_none = Malloc(char, 6);
	opt_p =  Malloc(char, 3);
	opt_ignh = Malloc(char, 7);
	opt_c = Malloc(char, 3);
	opt_rerun  = Malloc(char, 10);
	opt_e = Malloc(char, 3);
	opt_g = Malloc(char, 3);
	opt_num = Malloc(char, 5);
	g_sas_values = Malloc(double, MAX_VALUE_G_SAS);
	
	strcpy(opt_f, "-f");	
	strcpy(opt_o, "-o");
	strcpy(opt_ff, "-ff");
	strcpy(opt_water, "-water");
	strcpy(opt_none, "none");
	strcpy(opt_p, "-p");
	strcpy(opt_ignh, "-ignh");
	strcpy(opt_c, "-c");
	strcpy(opt_f, "-f");
	strcpy(opt_s, "-s");
	strcpy(opt_o, "-o");
	strcpy(f_step0, "step0*.pdb");
	strcpy(prot_gro, "prot.gro");
	strcpy(prot_top, "prot.top");
	strcpy(prot_tpr, "prot.tpr");
	strcpy(prot_trr, "prot.trr");
	strcpy(prot_log, "prot.log");
	strcpy(confout_gro, "confout.gro");
	strcpy(posre_itp, "posre.itp");
	strcpy(mdout_mdp, "mdout.mdp");
	strcpy(file_energy_computed_ener_edr, "file_energy_computed.ener.edr");
	strcpy(prot_sys_trr, "prot_sys.trr");
	strcpy(prot_sys_tpr, "prot_sys.tpr");
	strcpy(prot_sys_top, "prot_sys.top");
	strcpy(prot_sys_gro, "prot_sys.gro");	
	strcpy(energy_xvg, "energy.xvg");
	strcpy(traj_xtc, "traj.xtc");
	strcpy(xvg_1, "*.xvg_1*");
	strcpy(opt_rerun, "-rerun");
	strcpy(opt_e, "-e");	
	strcpy(opt_g, "-g");
	strcpy(opt_num, "-num");
	initialize_g_sas_values();

}

/** Finish GROMACS execution 
* This function is to be called when the Gromacs funcions below will no longer
* be used in order to release memory used by internal variables 
*/
_2PG_CARTESIAN_EXPORT
void finish_gromacs_execution(){

	free(command);
	free(program);
	free(filenm1);
	free(filenm2);
	free(filenm3);
	free(filenm4);
	free(filenm5);
	free(f_step0);
	free(prot_gro);
	free(prot_top);
	free(prot_tpr);
	free(prot_trr);
	free(prot_log);
	free(confout_gro);
	free(posre_itp);
	free(mdout_mdp);
	free(file_energy_computed_ener_edr);
	free(prot_sys_trr);
	free(prot_sys_tpr);
	free(prot_sys_top);
	free(prot_sys_gro);
	free(energy_xvg);
	free(traj_xtc);
	free(opt_f);
	free(opt_s);
	free(opt_o);
	free(opt_ff);
	free(opt_water);
	free(opt_none);
	free(opt_p);
	free(opt_ignh);
	free(opt_c);
	free(xvg_1);
	free(opt_rerun);
	free(opt_e);
	free(opt_g);
	free(opt_num);
	free(g_sas_values);
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

/** Calls toremove files used in simulation.
* Simulation means an execution to calculate the objectives by GROMACS.
*/
void clean_gromacs_simulation(const char *path_local_execute){
	delete_file(path_local_execute, f_step0);
	delete_file(path_local_execute, prot_gro);
	delete_file(path_local_execute, prot_top);
	delete_file(path_local_execute, prot_tpr);
	delete_file(path_local_execute, prot_trr);
	delete_file(path_local_execute, prot_log);	
	delete_file(path_local_execute, confout_gro);
	delete_file(path_local_execute, posre_itp);
	delete_file(path_local_execute, mdout_mdp);
	delete_file(path_local_execute, file_energy_computed_ener_edr);
	delete_file(path_local_execute, prot_sys_trr);
	delete_file(path_local_execute, prot_sys_tpr);
	delete_file(path_local_execute, prot_sys_top);
	delete_file(path_local_execute, prot_sys_gro);
	delete_file(path_local_execute, energy_xvg);
	delete_file(path_local_execute, traj_xtc);
	delete_file(path_local_execute, xvg_1);

}

/** Obtaing the value of objective in Double representation
* This function was created to be a pattern for getting
* the value of objective when it will be in Double
* instead of char.
*/
double get_objective_value(const char *value){
	double aux = MAX_ENERGY;
	aux = str2double(value);
	if ( (aux ==  -1.000000) || (isnan(aux) == btrue) ){
		aux = MAX_ENERGY;
	}
	return aux;
}

/* Sets the value of objective in solution 
* sol represents the solution (individual). The value of objective will be 
*     set in obj_values field.
* value the calculated value that will be set in solution
* obj represents the index of objective
*/
static void set_objective_from_gromacs_in_solution(solution_t *sol, char *value,
		const int *obj){
	sol->obj_values[*obj] = get_objective_value(value);
}


/** runs g_gyrate program and evaluates the radius of gyration of the protein
 * based on alpha-carbons only
 */
void compute_gyrate(solution_t *sol, const int *fit, const char *local_execute,
		const char *path_gromacs_programs, const char *pdbfile,
		const option_g_energy *opt_fitness, const char *computed_g_gyrate_value_file){	
	
	char *last_line, *line_splited;
	char *value;

	char *g_gyrate_args[8];

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
	strcat(filenm2, prot_tpr); /* local_execute/prot.tpr */
	g_gyrate_args[4] = filenm2;
	//xvg
	g_gyrate_args[5] = opt_o;
	strcpy(filenm3, local_execute);
	strcat(filenm3, computed_g_gyrate_value_file);	
	g_gyrate_args[6] = filenm3;
	//final parameter
	g_gyrate_args[7] = NULL;

	if (!run_program_after_pipe("C-alpha", program, g_gyrate_args))
		fatal_error("Failed to run g_gyrate at compute_gyrate function\n");

	if (check_exists_file(computed_g_gyrate_value_file) == btrue){
		value = Malloc(char,MAX_VALUE);
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
	    delete_file(local_execute, computed_g_gyrate_value_file);
	}else{
		sol->obj_values[*fit] = MAX_ENERGY; //MAX energy value
	}
}

static void build_pdb_file_name(char *pdb_file_name, const char *aux_name,
		const char *__restrict prefix){
	strcpy(pdb_file_name, prefix);
	strcat(pdb_file_name,aux_name);
	strcat(pdb_file_name,".pdb");
}

/** Creates a tpr file
*/
void build_tpr_file(const char *pdbfile, const char *local_execute,
		const char *path_gromacs_programs, const char *force_field, const char *mdp_file){

	const protein_t *protein_aux;	
	
	char *pdbfile_aux;
	char *force_field_aux;
	char *mdp_file_aux;

	char *pdb2gmx_args[13];
	char *grompp_args[11];

	pdbfile_aux = Malloc(char, MAX_FILE_NAME);
	force_field_aux = Malloc(char, MAX_FORCE_FIELD_NAME);
	mdp_file_aux = Malloc(char, MAX_FILE_NAME);

	/* pdb2gmx */
	strcpy(program, path_gromacs_programs);
	strcat(program, "pdb2gmx");
	pdb2gmx_args[0] = program;
	//pdb
	strcpy(pdbfile_aux, pdbfile);
	pdb2gmx_args[1] = opt_f;
	strcpy(filenm1, local_execute);	
	strcat(filenm1, pdbfile_aux);
	pdb2gmx_args[2] = filenm1;
	//gro
	pdb2gmx_args[3] = opt_o;
	strcpy(filenm2, local_execute);
	strcat(filenm2, "prot.gro");
	pdb2gmx_args[4] = filenm2;
	//top
	pdb2gmx_args[5] = opt_p;
	strcpy(filenm3, local_execute);
	strcat(filenm3, "prot.top");
	pdb2gmx_args[6] = filenm3;
	//force field
	pdb2gmx_args[7] = opt_ff;
	strcpy(force_field_aux, force_field);
	pdb2gmx_args[8] = force_field_aux;
	//water
	pdb2gmx_args[9] = opt_water;
	pdb2gmx_args[10] = opt_none;
	//Hydrogen
	pdb2gmx_args[11] = opt_ignh;

	pdb2gmx_args[12] = NULL;	

	if (!run_program(program, pdb2gmx_args)){
		fatal_error("Failed to run pdb2gmx at build_tpr_file function \n");
	}

	/* grompp */
	strcpy(program, path_gromacs_programs);
	strcat(program, "grompp");
	grompp_args[0] = program;
	//mdp
	grompp_args[1] = opt_f;
	strcpy(mdp_file_aux, mdp_file);
	strcpy(filenm1, local_execute);
	strcat(filenm1, mdp_file_aux);
	grompp_args[2] = filenm1;
	//top
	grompp_args[3] = opt_p;
	strcpy(filenm2, local_execute);
	strcat(filenm2, "prot.top");
	grompp_args[4] = filenm2;
	//tpŕ
	grompp_args[5] = opt_o;
	strcpy(filenm3, local_execute);
	strcat(filenm3, "prot.tpr");
	grompp_args[6] = filenm3;
	//gro
	grompp_args[7] = opt_c;
	strcpy(filenm4, local_execute);
	strcat(filenm4, "prot.gro");
	grompp_args[8] = filenm4;

	grompp_args[9] = NULL;
	grompp_args[10] = NULL;

	if (!run_program(program, grompp_args))
		fatal_error("Failed to run grompp at build_tpr_file function \n");

	free(pdbfile_aux);
	free(force_field_aux);
	free(mdp_file_aux);
}

/** Calls mdrun program to calculate energies
* Energies will be used as objectivies
* It execution was based on http://www.gromacs.org/Documentation/How-tos/Single-Point_Energy
*/
void call_mdrun2energy(const char *pdbfile, const char *local_execute,
		const char *path_gromacs_programs){
	char *mdrun_args[13];
	/* mdrun */
	strcpy(program, path_gromacs_programs);
	strcat(program, "mdrun");
	mdrun_args[0] = program;
	//tpŕ
	mdrun_args[1] = opt_s;		
	strcpy(filenm1, local_execute);
	strcat(filenm1, prot_tpr);
	mdrun_args[2] = filenm1;
	//pdb
	mdrun_args[3] = opt_rerun;
	strcpy(filenm2, local_execute);
	strcat(filenm2, pdbfile);
	mdrun_args[4] = filenm2;
	//trr
	mdrun_args[5] = opt_o;
	strcpy(filenm3, local_execute);
	strcat(filenm3, prot_trr);
	mdrun_args[6] = filenm3;
	//energy file
	mdrun_args[7] = opt_e;
	strcpy(filenm4, local_execute);
	strcat(filenm4, file_energy_computed_ener_edr);
	mdrun_args[8] = filenm4;
	//log file
	mdrun_args[9] = opt_g;
	strcpy(filenm5, local_execute);
	strcat(filenm5, prot_log);
	mdrun_args[10] = filenm5;

	mdrun_args[11] = NULL;
	mdrun_args[12] = NULL;

	if (!run_program(program, mdrun_args))
		fatal_error("Failed to run mdrun at call_mdrun2energy function \n");
}

/** Checks the opt_objective is one energy objective
*/
boolean_t is_energy_objective(const option_g_energy *opt_objective){

	if ( (*opt_objective == gmx_elel_ener) ||
		 (*opt_objective == gmx_edw_ener) ||
		 (*opt_objective == gmx_potential_ener)	||
		 (*opt_objective == gmx_GBSA_Solvatation) ){
		 	return btrue;
	}
	return bfalse;
}	

/** Calls the g_energy program
*/
void call_g_energy(const char *local_execute, const char *path_gromacs_programs, 
	const char *opt_energy ){

	char *g_energy_args[7];

	/*g_energy*/	
	strcpy(program, path_gromacs_programs);
	strcat(program, "g_energy");
	g_energy_args[0] = program;
	//Energy file
	g_energy_args[1] = opt_f;
	strcpy(filenm1, local_execute);
	strcat(filenm1, file_energy_computed_ener_edr);
	g_energy_args[2] = filenm1;
	//xvg
	g_energy_args[3] = opt_o;
	strcpy(filenm2, local_execute);
	strcat(filenm2, energy_xvg);	
	g_energy_args[4] = filenm2;

	g_energy_args[5] = NULL;
	g_energy_args[6] = NULL;
	
	if (!run_program_after_pipe(opt_energy, program, g_energy_args))
		fatal_error("Failed to run g_energy at compute_energy_minimum function\n");


}

/** Calculates the energies, except GBSA Solvatation
*/
void compute_energy(solution_t *sol, const int *obj, const char *local_execute,
		const char *path_gromacs_programs, const char *opt_energy ){
	char *last_line, *line_splited;
	char *value;
	
	if (check_exists_file(file_energy_computed_ener_edr) == btrue){
		//Call g_energy
		call_g_energy(local_execute, path_gromacs_programs, opt_energy);	
		if (check_exists_file(energy_xvg) == btrue){
			value = Malloc(char,MAX_VALUE);
			//get the last line of xvg file
			last_line = get_last_line(energy_xvg);
			// Split last line by space and obtaing the first value
	 		line_splited = strtok (last_line," ");
	 		// Obtaining the second value. It will be set in solution
	 		line_splited = strtok(NULL, " ");
	 		strcpy(value, line_splited);
	 		//Looking the end of line_splited
		  	while (line_splited != NULL){	    	
	    		line_splited = strtok(NULL, " ");
	  		}
			//Set the value of energy option in solution
		    set_objective_from_gromacs_in_solution(sol,value, obj);	
		    free(line_splited);
		    free(last_line);
		    free(value);
		    delete_file(local_execute, energy_xvg);
		}else{
			sol->obj_values[*obj] = MAX_ENERGY; //MAX energy value
		}
	}else{
		sol->obj_values[*obj] = MAX_ENERGY; //MAX energy value
	}

}

/** Calculates GBSA Solvatation energy
* It value is composed by two values: gmx_GB_Polarization and gmx_Nonpolar_Sol
* Therefore g_energy program is called twice. The value of gmx_GB_Polarization
* is stored at G_pol. the value of gmx_Nonpolar_Sol is stored at G_np. Therefore,
* the value of GBSA is stored at G_solv. 
* G_solv variable is setted in solution.
* Important: When error is found either gmx_GB_Polarization or gmx_Nonpolar_Sol 
* solution will set MAX value
*/
void compute_energy_GBSA(solution_t *sol, const int *obj, const char *local_execute,
		const char *path_gromacs_programs, const char *opt_energy ){
	char *last_line, *line_splited;
	char *value;
	double G_solv, G_pol, G_np;	

	//Call g_energy - gmx_GB_Polarization
	call_g_energy(local_execute, path_gromacs_programs, 
		option_g_energy_program[gmx_GB_Polarization].option_name);
	if (check_exists_file(energy_xvg) == btrue){
		value = Malloc(char,MAX_VALUE);
		//get the last line of xvg file
		last_line = get_last_line(energy_xvg);
		// Split last line by space and obtaing the first value
 		line_splited = strtok (last_line," ");
 		// Obtaining the second value. It will be set in solution
 		line_splited = strtok(NULL, " ");
 		strcpy(value, line_splited);
 		//Looking the end of line_splited
	  	while (line_splited != NULL){	    	
    		line_splited = strtok(NULL, " ");
  		}
		//Set the value of energy option in solution
	    G_pol = get_objective_value(value);
	    free(line_splited);
	    free(last_line);
	    free(value);
	    delete_file(local_execute, energy_xvg);
		//Call g_energy - gmx_Nonpolar_Sol
		call_g_energy(local_execute, path_gromacs_programs, 
			option_g_energy_program[gmx_Nonpolar_Sol].option_name);
		if (check_exists_file(energy_xvg) == btrue){
			value = Malloc(char,MAX_VALUE);
			//get the last line of xvg file
			last_line = get_last_line(energy_xvg);
			// Split last line by space and obtaing the first value
 			line_splited = strtok (last_line," ");
 			// Obtaining the second value. It will be set in solution
 			line_splited = strtok(NULL, " ");
 			strcpy(value, line_splited);
 			//Looking the end of line_splited
	  		while (line_splited != NULL){	    	
    			line_splited = strtok(NULL, " ");
  			}
			//Set the value of energy option in solution
	    	G_np = get_objective_value(value);
			//Getting the value of GBSA Solvatation
			G_solv = G_pol + G_np;
			sol->obj_values[*obj] = G_solv;
	    	free(line_splited);
	    	free(last_line);
	    	free(value);	    
	    	delete_file(local_execute, energy_xvg);
		}else{
			sol->obj_values[*obj] = MAX_ENERGY;
		}
	}else{
		sol->obj_values[*obj] = MAX_ENERGY;
	}
}

/** Checks the opt_objective is one sas objective
*/
boolean_t is_sas_objective(const option_g_energy *opt_objective){

	if ( (*opt_objective == gmx_hydrophobic) ||
		 (*opt_objective == gmx_hydrophilic) ||
		 (*opt_objective == gmx_total_area) ){
		 	return btrue;
	}
	return bfalse;
}

/** runs g_sas program to calculate hydrophobic, hydrophilic and total 
* solvent accessible surface area. 
* See Eisenhaber F, Lijnzaad P, Argos P, Sander C, & Scharf M (1995) J. Comput. Chem. 16, 273-284.
 */
void call_g_sas(const char *local_execute, const char *path_gromacs_programs, 
	const char *pdbfile, const char *computed_g_sas_value_file){
	char *g_sas_args[9];
	char *last_line, *line_splited;
	char *value;

	/* g_sas */
	strcpy(program, path_gromacs_programs);
	strcat(program, "g_sas");
	g_sas_args[0] = program;
	//gro instead of PDB
	g_sas_args[1] = opt_f;
	strcpy(filenm1, local_execute);
	strcat(filenm1, prot_gro);
	g_sas_args[2] = filenm1;
	//tpr
	g_sas_args[3] = opt_s;
	strcpy(filenm2, local_execute);
	strcat(filenm2, prot_tpr);
	g_sas_args[4] = filenm2;
	//xvg
	g_sas_args[5] = opt_o;
	strcpy(filenm3, local_execute);	
	strcat(filenm3, computed_g_sas_value_file);
	g_sas_args[6] = filenm3;

	g_sas_args[7] = NULL;
	g_sas_args[8] = NULL;

	if (!run_program_after_pipe("Protein Protein", program, g_sas_args)){
		fatal_error("Failed to run g_sas program at call_g_sas function\n");
	}		
	if (check_exists_file(computed_g_sas_value_file)){
		value = Malloc(char,MAX_VALUE);
		//get the last line of xvg file
		last_line = get_last_line(computed_g_sas_value_file);
		// Split last line by space and obtaing the first value
 		line_splited = strtok (last_line," ");
 		//Set the value for hydrophobic objective		
 		line_splited = strtok(NULL, " ");
 		strcpy(value, line_splited);
 		g_sas_values[0] = get_objective_value(value);
 		//Set the value for hydrophilic objective
 		line_splited = strtok(NULL, " ");
 		strcpy(value, line_splited);
 		g_sas_values[1] = get_objective_value(value);
 		g_sas_values[1] = g_sas_values[1] * (-1);
 		//Set the value for total area objective
 		line_splited = strtok(NULL, " ");
 		strcpy(value, line_splited);
 		g_sas_values[2] = get_objective_value(value); 		
 		//Looking the end of line_splited
	  	while (line_splited != NULL){	    	
    		line_splited = strtok(NULL, " ");
  		}	    	
	    free(line_splited);
	    free(last_line);
	    free(value);
		delete_file(local_execute, computed_g_sas_value_file);
	}else{
		//Hydrophobic
		g_sas_values[0] = MAX_ENERGY;
		//Hydrophilic
		g_sas_values[1] = MIN_ENERGY;
		//Total Area 
		g_sas_values[3] = MAX_ENERGY;
	}
}

/** Runs g_hbond program to compute the number of hydrogen bonds of the protein 
* The Group is Protein Protein
*/
void compute_hbond(solution_t *sol, const int *fit, const char *local_execute,
		const char *path_gromacs_programs, const char *pdbfile,
		const option_g_energy *opt_fitness, const char *computed_hbond_value_file){	
	
	char *last_line, *line_splited;
	char *value;

	char *g_hbond_args[9];

	/* g_hbond */
	strcpy(program, path_gromacs_programs);
	strcat(program, "g_hbond");
	g_hbond_args[0] = program;
	//gro instead of pdb
	g_hbond_args[1] = opt_f;
	strcpy(filenm1, local_execute);
	strcat(filenm1, prot_gro);
	g_hbond_args[2] = filenm1;
	//tpr
	g_hbond_args[3] = opt_s;
	strcpy(filenm2, local_execute);
	strcat(filenm2, prot_tpr);
	g_hbond_args[4] = filenm2;
	//xvg
	g_hbond_args[5] = opt_num;
	strcpy(filenm3, local_execute);
	strcat(filenm3, computed_hbond_value_file);
	g_hbond_args[6] = filenm3;

	g_hbond_args[7] = NULL;
	g_hbond_args[8] = NULL;

	if (!run_program_after_pipe("Protein Protein", program, g_hbond_args))
		fatal_error("Failed to run g_hbond at compute_hbond function\n");

	if (check_exists_file(computed_hbond_value_file) == btrue){
		value = Malloc(char,MAX_VALUE);
		//get the last line of xvg file
		last_line = get_last_line(computed_hbond_value_file);
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
	    //It will be maximized that means minimize its opposite value
	    sol->obj_values[*fit] = sol->obj_values[*fit] * (-1);
	    free(line_splited);
	    free(last_line);
	    free(value);
	    delete_file(local_execute, computed_hbond_value_file);
	}else{
		sol->obj_values[*fit] = MIN_ENERGY;
	}
}

/** Runs g_hbond program to compute the number of hydrogen bonds of the protein 
* The Group is "MainChain+H MainChain+H"
*/
void compute_hbond_main(solution_t *sol, const int *fit, const char *local_execute,
		const char *path_gromacs_programs, const char *pdbfile,
		const option_g_energy *opt_fitness, const char *computed_hbond_value_file){	
	
	char *last_line, *line_splited;
	char *value;

	char *g_hbond_args[9];

	/* g_hbond */
	strcpy(program, path_gromacs_programs);
	strcat(program, "g_hbond");
	g_hbond_args[0] = program;
	//gro instead of pdb
	g_hbond_args[1] = opt_f;
	strcpy(filenm1, local_execute);
	strcat(filenm1, prot_gro);
	g_hbond_args[2] = filenm1;
	//tpr
	g_hbond_args[3] = opt_s;
	strcpy(filenm2, local_execute);
	strcat(filenm2, prot_tpr);
	g_hbond_args[4] = filenm2;
	//xvg
	g_hbond_args[5] = opt_num;
	strcpy(filenm3, local_execute);
	strcat(filenm3, computed_hbond_value_file);
	g_hbond_args[6] = filenm3;

	g_hbond_args[7] = NULL;
	g_hbond_args[8] = NULL;

	if (!run_program_after_pipe("MainChain+H MainChain+H", program, g_hbond_args))
		fatal_error("Failed to run g_hbond at compute_hbond_main function\n");

	if (check_exists_file(computed_hbond_value_file) == btrue){
		value = Malloc(char,MAX_VALUE);
		//get the last line of xvg file
		last_line = get_last_line(computed_hbond_value_file);
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
	    //It will be maximized that means minimize its opposite value
	    sol->obj_values[*fit] = sol->obj_values[*fit] * (-1);
	    free(line_splited);
	    free(last_line);
	    free(value);
	    delete_file(local_execute, computed_hbond_value_file);
	}else{
		sol->obj_values[*fit] = MIN_ENERGY;
	}
}

/** Initialize values for next solution
* This function is used to guarantee that the values for next
* solution will be initialized
*/
void initialize_values_for_next_solution(){
	initialize_g_sas_values();
}

/** Calculates the objectives by GROMACS
*/
_2PG_CARTESIAN_EXPORT
void get_gromacs_objectives(solution_t *solutions, const input_parameters_t *in_para){
	const protein_t *population_aux;
	char pdbfile_aux[30];
	char aux [20];
	char aux_ind[5];
	char msg [50];	
	option_g_energy *opt_objective;
	int has_energy_objective;
	int has_sas_objective;

	has_energy_objective = 0;
	has_sas_objective = 0;
	//getting the gromacs objectives from input parameters
	opt_objective = Malloc(option_g_energy,in_para->number_fitness);
	for (int ob = 0; ob < in_para->number_fitness; ob++){
		opt_objective[ob] = get_option_g_energy_t_from_type_fitness_energy(&in_para->fitness_energies[ob]);
		if ( is_energy_objective(&opt_objective[ob]) == btrue){
			has_energy_objective = has_energy_objective + 1;
		}else  if ( is_sas_objective(&opt_objective[ob]) == btrue){
			has_sas_objective = has_sas_objective + 1;
		}
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
   	    //Building Generic tpr file. It will be used in all GROMACS execution   	    
    	build_tpr_file(pdbfile_aux, in_para->path_local_execute, in_para->path_gromacs_programs, 
	        	in_para->force_field, in_para->mdp_file);
    	/*Check if call mdrun to calculate the energies. 
    	 * When has_energy_objective is greater than 0, it means one energy objective exists at least 
    	*/
    	if (has_energy_objective > 0){
    		call_mdrun2energy(pdbfile_aux, in_para->path_local_execute, 
	    		in_para->path_gromacs_programs);
    	}
    	if (has_sas_objective > 0){
			call_g_sas(in_para->path_local_execute, in_para->path_gromacs_programs, pdbfile_aux, 
				in_para->computed_areas_g_sas_file);
    	}
   	    // obtaing the values of objectivies
   	    for (int ob = 0; ob < in_para->number_fitness; ob++){
	   	    if (opt_objective[ob] == gmx_gyrate) {
   		    	compute_gyrate(&solutions[ind], &ob, in_para->path_local_execute,
   	    				in_para->path_gromacs_programs, pdbfile_aux, opt_objective,
   	    				in_para->computed_radius_g_gyrate_file);
   	    	}else if ( (opt_objective[ob] == gmx_potential_ener) ||
   	    	 			(opt_objective[ob] == gmx_elel_ener) ||
   	    	 			(opt_objective[ob] == gmx_edw_ener)){   	    		
   	    		compute_energy(&solutions[ind], &ob, in_para->path_local_execute,
   	    				in_para->path_gromacs_programs,
   	    				option_g_energy_program[opt_objective[ob]].value_opt );
   	    	}else if ( opt_objective[ob] == gmx_GBSA_Solvatation){   	    		
   	    		compute_energy_GBSA(&solutions[ind], &ob, in_para->path_local_execute,
   	    				in_para->path_gromacs_programs,
   	    				option_g_energy_program[opt_objective[ob]].value_opt );
   	    	}else if ( opt_objective[ob] == gmx_hydrophobic){
   	    		solutions[ind].obj_values[ob] = g_sas_values[0];
   	    	}else if ( opt_objective[ob] == gmx_hydrophilic){
   	    		solutions[ind].obj_values[ob] = g_sas_values[1];
   	    	}else if ( opt_objective[ob] == gmx_total_area){
   	    		solutions[ind].obj_values[ob] = g_sas_values[2];
   	    	}else if ( opt_objective[ob] == gmx_hbond){
   		    	compute_hbond(&solutions[ind], &ob, in_para->path_local_execute,
   	    				in_para->path_gromacs_programs, pdbfile_aux, opt_objective,
   	    				in_para->computed_g_hbond_file);   	    		
   	    	}else if ( opt_objective[ob] == gmx_hbond_main){
   		    	compute_hbond_main(&solutions[ind], &ob, in_para->path_local_execute,
   	    				in_para->path_gromacs_programs, pdbfile_aux, opt_objective,
   	    				in_para->computed_g_hbond_file);   	    		
   	    	}
   	    }
   	    initialize_values_for_next_solution();
   	    clean_gromacs_simulation(in_para->path_local_execute);
	}

	free(opt_objective);
}


/** Calculates the objectives of one solution by GROMACS
*/
_2PG_CARTESIAN_EXPORT
void get_gromacs_objectives_of_solution(solution_t *solution, 
	const input_parameters_t *in_para, const int *ind){
	const protein_t *population_aux;
	char pdbfile_aux[30];
	char aux [20];
	char aux_ind[5];
	char msg [50];	
	option_g_energy *opt_objective;
	int has_energy_objective;
	int has_sas_objective;

	has_energy_objective = 0;
	has_sas_objective = 0;
	//getting the gromacs objectives from input parameters
	opt_objective = Malloc(option_g_energy,in_para->number_fitness);
	for (int ob = 0; ob < in_para->number_fitness; ob++){
		opt_objective[ob] = get_option_g_energy_t_from_type_fitness_energy(&in_para->fitness_energies[ob]);
		if ( is_energy_objective(&opt_objective[ob]) == btrue){
			has_energy_objective = has_energy_objective + 1;
		}else  if ( is_sas_objective(&opt_objective[ob]) == btrue){
			has_sas_objective = has_sas_objective + 1;
		}
	}
	//getting the protein from solution
	population_aux = (protein_t*) solution->representation;
	// saving pdb file of protein
	int2str(aux_ind,ind);
	build_pdb_file_name(pdbfile_aux, aux_ind, PREFIX_PDB_FILE_NAME_EA);
	save_pdb_file(in_para->path_local_execute, pdbfile_aux, 
	    	&population_aux->p_topol->numatom, population_aux->p_atoms, NULL);
	    //Building Generic tpr file. It will be used in all GROMACS execution   	    
	build_tpr_file(pdbfile_aux, in_para->path_local_execute, in_para->path_gromacs_programs, 
        	in_para->force_field, in_para->mdp_file);
	/*Check if call mdrun to calculate the energies. 
	 * When has_energy_objective is greater than 0, it means one energy objective exists at least 
	*/
	if (has_energy_objective > 0){
		call_mdrun2energy(pdbfile_aux, in_para->path_local_execute, 
    		in_para->path_gromacs_programs);
	}
	if (has_sas_objective > 0){
		call_g_sas(in_para->path_local_execute, in_para->path_gromacs_programs, pdbfile_aux, 
			in_para->computed_areas_g_sas_file);
	}
    // obtaing the values of objectivies
    for (int ob = 0; ob < in_para->number_fitness; ob++){
	    if (opt_objective[ob] == gmx_gyrate) {
	    	compute_gyrate(solution, &ob, in_para->path_local_execute,
    				in_para->path_gromacs_programs, pdbfile_aux, opt_objective,
    				in_para->computed_radius_g_gyrate_file);
    	}else if ( (opt_objective[ob] == gmx_potential_ener) ||
    	 			(opt_objective[ob] == gmx_elel_ener) ||
    	 			(opt_objective[ob] == gmx_edw_ener)){   	    		
    		compute_energy(solution, &ob, in_para->path_local_execute,
    				in_para->path_gromacs_programs,
    				option_g_energy_program[opt_objective[ob]].value_opt );
    	}else if ( opt_objective[ob] == gmx_GBSA_Solvatation){   	    		
    		compute_energy_GBSA(solution, &ob, in_para->path_local_execute,
    				in_para->path_gromacs_programs,
    				option_g_energy_program[opt_objective[ob]].value_opt );
    	}else if ( opt_objective[ob] == gmx_hydrophobic){
    		solution->obj_values[ob] = g_sas_values[0];
    	}else if ( opt_objective[ob] == gmx_hydrophilic){
    		solution->obj_values[ob] = g_sas_values[1];
    	}else if ( opt_objective[ob] == gmx_total_area){
    		solution->obj_values[ob] = g_sas_values[2];
    	}else if ( opt_objective[ob] == gmx_hbond){
	    	compute_hbond(solution, &ob, in_para->path_local_execute,
    				in_para->path_gromacs_programs, pdbfile_aux, opt_objective,
    				in_para->computed_g_hbond_file);   	    		
    	}else if ( opt_objective[ob] == gmx_hbond_main){
	    	compute_hbond_main(solution, &ob, in_para->path_local_execute,
    				in_para->path_gromacs_programs, pdbfile_aux, opt_objective,
    				in_para->computed_g_hbond_file);   	    		
    	}
    }
    initialize_values_for_next_solution();
    clean_gromacs_simulation(in_para->path_local_execute);
	
	free(opt_objective);
}


/** Run of pdb2gmx for pattern of atom names
* Output is the same name of pdbfile. 
* WARNING: pdbfile will be written
*/
void call_pdb2gmx_for_pattern_atom_names(const char *pdbfile, const char *local_execute,
		const char *path_gromacs_programs, const char *force_field){

	const protein_t *protein_aux;	
	
	char *pdbfile_aux;
	char *force_field_aux;

	char *pdb2gmx_args[13];
	char *grompp_args[11];

	pdbfile_aux = Malloc(char, MAX_FILE_NAME);
	force_field_aux = Malloc(char, MAX_FORCE_FIELD_NAME);

	/* pdb2gmx */
	strcpy(program, path_gromacs_programs);
	strcat(program, "pdb2gmx");
	pdb2gmx_args[0] = program;
	//pdb
	strcpy(pdbfile_aux, pdbfile);
	pdb2gmx_args[1] = opt_f;
	strcpy(filenm1, local_execute);	
	strcat(filenm1, pdbfile_aux);
	pdb2gmx_args[2] = filenm1;
	//gro
	pdb2gmx_args[3] = opt_o;
	strcpy(filenm2, local_execute);
	strcat(filenm2, pdbfile_aux);
	pdb2gmx_args[4] = filenm2;
	//top
	pdb2gmx_args[5] = opt_p;
	strcpy(filenm3, local_execute);
	strcat(filenm3, "prot.top");
	pdb2gmx_args[6] = filenm3;
	//force field
	pdb2gmx_args[7] = opt_ff;
	strcpy(force_field_aux, force_field);
	pdb2gmx_args[8] = force_field_aux;
	//water
	pdb2gmx_args[9] = opt_water;
	pdb2gmx_args[10] = opt_none;
	//Hydrogen
	pdb2gmx_args[11] = opt_ignh;

	pdb2gmx_args[12] = NULL;	

	if (!run_program(program, pdb2gmx_args)){
		fatal_error("Failed to run pdb2gmx at call_pdb2gmx_for_pattern_atom_names function \n");
	}

	free(pdbfile_aux);
	free(force_field_aux);
}
