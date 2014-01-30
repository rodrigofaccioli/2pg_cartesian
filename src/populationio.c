#include<stdio.h>
#include<string.h>

//#include "protein.h"
#include "populationio.h"
#include "futil.h"
#include "math_owner.h"
#include <stdlib.h>

static void start_ind(FILE *pop_file){
	fprintf (pop_file,"##\n");
}

static void finish_ind(FILE *pop_file){
	fprintf (pop_file,"$$\n");
}

static void write_individual(FILE *pop_file,protein **pop,
		const int *ind){
	write_protein(pop_file,pop[*ind]);
}

static void write_protein(FILE *pop_file, const protein *prot){
	for (int r = 0; r < prot->nr_residues; r++){
		fprintf(pop_file,"%f %f %f %d %d",radians2degree(&prot->residuo[r].phi) ,
				radians2degree(&prot->residuo[r].psi), radians2degree(&prot->residuo[r].omega),
				prot->residuo[r].number_late, prot->residuo[r].pos_late);
		for (int s = 0; s < prot->residuo[r].number_late;s++){
			fprintf(pop_file," %f", radians2degree(&prot->residuo[r].late[s]));
		}
		fprintf (pop_file," \n");
	}
}
void _save_population_file(const char *path, const char *file_name,
		protein ** pop, const int *pop_size){
	write_initial_population_file(path, file_name,pop, pop_size);
}

void _save_protein_path_file(const char *path_file_name, const protein *prot){
	/*Save a protein. Therefore, the individual number is 1 */
	FILE *prot_file = open_file(path_file_name,fWRITE);
	int ind = 0;

	start_ind(prot_file);
	write_protein(prot_file,&prot[ind]);
	finish_ind(prot_file);

	fclose(prot_file);

}

void write_initial_population_file(const char *path, const char *file_name,
		protein ** pop, const int *pop_size){
        char *fname = path_join_file(path,file_name);
	FILE *pop_file = open_file(fname,fWRITE);

	for (int ind = 0; ind < *pop_size; ind++){
		start_ind(pop_file);
		write_individual(pop_file,pop,&ind);
		finish_ind(pop_file);
	}
	fclose(pop_file);
	free(fname);
}

