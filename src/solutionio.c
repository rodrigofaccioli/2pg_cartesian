#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "solution.h"
#include "defines.h"
#include "solutionio.h"
#include "futil.h"
#include "string_owner.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

#define MAX_LINE_SOLUTION 200

static void write_header_generation(FILE *fit_file, const int *ger){
	fprintf (fit_file,"#Generation %d \n", *ger);
	fprintf (fit_file,"#Solution Objective \n");

}

static void write_objectives_values(FILE *fit_file, const int *fit, 
	const solution_t *solutions, const int *pop_size){
	int i;
	for (i=0; i < *pop_size; i++){
		fprintf (fit_file,"%d %f \n", i+1, get_specific_objective(&solutions[i],fit));
	}
}

static void write_oposite_objective_values(FILE *fit_file, const int *fit, 
	const solution_t *solutions, const int *pop_size){
	int i;
	for (i=0; i < *pop_size; i++){
		fprintf (fit_file,"%d %f \n", i+1, get_oposite_specific_objective(&solutions[i],fit));
	}
}

void save_solution_file(const char *path, const char *file_name,
		const int *fit, const solution_t *solutions, const int *pop_size,
		const int *tag, const input_parameters_t *in_para){

    char *fname = path_join_file(path,file_name);
	FILE *fit_file = open_file(fname,fWRITE);
	write_header_generation(fit_file, tag);
	/* objective that must be maximized. When they are obtained these values
	 * are multiplied by -1 because 2PG works with these opposite values.
	 * However, when they will be stored, they must be written in original
	 * value */
	if ( (in_para->fitness_energies[*fit] == fit_hbond) ||
		(in_para->fitness_energies[*fit] == fit_hydrophilic)||
		(in_para->fitness_energies[*fit] == fit_hbond_main)	||
		(in_para->fitness_energies[*fit] == fit_stride_total) ||
		(in_para->fitness_energies[*fit] == fit_stride_helix)	||
		(in_para->fitness_energies[*fit] == fit_stride_beta) ) {
		write_oposite_objective_values(fit_file,fit,solutions,pop_size);
	}else{
		write_objectives_values(fit_file,fit,solutions,pop_size);
	}
	fclose(fit_file);
	free(fname);
}


static void initialize_for_reading_file(int *num_solutions_r,
	int *numobj_r, const char *path_file_name){
	FILE * solution_file;	
	char *line;
	char *line_splited;
	int num_obj;
	int num_solutions;
	

	num_solutions = 0;
	num_obj = -1; //Because the first collum is index.
	line = Malloc(char, MAX_LINE_SOLUTION);
	solution_file = open_file(path_file_name, fREAD);

	//getting first line to obtain number of objectivies
	fgets(line,MAX_LINE_SOLUTION,solution_file);
 	line_splited = strtok (line,"\t");
  	while (line_splited != NULL){
    	num_obj = num_obj + 1;
    	line_splited = strtok(NULL, "\t");
  	}	
	//Counting how many lines the file has
	while ( fgets(line,MAX_LINE_SOLUTION,solution_file) != NULL){
		num_solutions = num_solutions + 1;
	}

	fclose(solution_file);
	free(line);

	//Setting return values
	*numobj_r = num_obj;
	*num_solutions_r = num_solutions;

}

/** Loading a solution file.

* Returns: an array of solution
* Impotant: This array of solution need	to be desallocated
*/
_2PG_CARTESIAN_EXPORT
solution_t * loading_file_solutions(int *num_solutions_r, 
	int *numobj_r, 	const char *path_file_name){
	FILE * solution_file;	
	char *line;
	char *line_splited;
	int sol;
	solution_t * solutions_aux;

	line = Malloc(char, MAX_LINE_SOLUTION);

	//Getting initial information
	initialize_for_reading_file(num_solutions_r, numobj_r, path_file_name);

	//Alocating Solution
	solutions_aux = allocate_solution(num_solutions_r, numobj_r);

	//Reading file and set values of objective
	sol = -1;
	solution_file = open_file(path_file_name, fREAD);
	//Removing first line that is collumn
	fgets(line,MAX_LINE_SOLUTION,solution_file);
	while ( fgets(line,MAX_LINE_SOLUTION,solution_file) != NULL){
		sol = sol + 1;
		//Removing index collumn
		line_splited = strtok (line,"\t");
		//Setting number of objectives
		solutions_aux[sol].num_obj = *numobj_r;
		for (int ob = 0; ob < *numobj_r; ob++){
			line_splited = strtok (NULL,"\t");
			trim(line_splited);
			solutions_aux[sol].obj_values[ob] = str2float(line_splited);//str2double(line_splited);
		}					
	}	
	fclose(solution_file);

	free(line);

	return solutions_aux;
}
