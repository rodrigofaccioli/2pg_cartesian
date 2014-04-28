#include <stdio.h>
#include <stdlib.h>

#include "solution.h"
#include "solutionio.h"
#include "futil.h"

static void write_header_generation(FILE *fit_file, const int *ger){
	fprintf (fit_file,";Generation %d \n", *ger);
	fprintf (fit_file,";Solution Objective \n");

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
