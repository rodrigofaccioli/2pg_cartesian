#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dominance.h"
#include "defines.h"
#include "messages.h"
#include "futil.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif


_2PG_CARTESIAN_EXPORT
dominance_t * allocate_dominance(const int *size){
	dominance_t *dominance;
	dominance = Malloc(dominance_t, *size);
	return dominance;
}

_2PG_CARTESIAN_EXPORT
void desallocate_dominance(dominance_t *dominance, const int *size){
	for (int d = 0; d < *size; d++){
		free(dominance[d].set_dominated);
	}	
	free(dominance);
	dominance = NULL;
}

_2PG_CARTESIAN_EXPORT
void show_dominance(const dominance_t *dominance, const int *size){
	for (int d = 0; d < *size; d++){
		printf("Index dominance %i\n", d+1);
		printf("How many solutions dominate me %i\n", 
			dominance[d].how_many_solutions_dominate_it);
		printf("How many solutions are dominated by me %i\n", 
			dominance[d].max_dominated);
		printf("Values of objective\n");
		for (int ob = 0; ob < dominance[d].sol->num_obj; ob++){
			printf("objective %i value is %f\n", ob, 
				dominance[d].sol->obj_values[ob]);
		}
		printf("-------------------- \n\n");
	}
}

_2PG_CARTESIAN_EXPORT
void save_dominance(const dominance_t *dominance, const int *size){
	FILE *d_file=NULL;
	char *file_name;
	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "dominance_output.txt");	
    
	d_file = open_file(file_name, fWRITE);

	for (int d = 0; d < *size; d++){
		fprintf(d_file, "Index dominance %i\n", d+1);
		fprintf(d_file, "How many solutions dominate me %i\n", 
			dominance[d].how_many_solutions_dominate_it);
		fprintf(d_file, "How many solutions are dominated by me %i\n", 
			dominance[d].max_dominated);
		fprintf(d_file, "Values of objective\n");
		for (int ob = 0; ob < dominance[d].sol->num_obj; ob++){
			fprintf(d_file, "objective %i value is %f\n", ob, 
				dominance[d].sol->obj_values[ob]);
		}
		fprintf(d_file,"-------------------- \n\n");
	}
	fclose(d_file);
	free(file_name);	
}


/** Applies the dominance concept in solutions
* dominance is a pointer that stores the application of dominance concept in solution
* solutions array of Solutions
* size number of solutions
*/
_2PG_CARTESIAN_EXPORT
void set_dominance(dominance_t *dominance, const solution_t *solutions, const int *size){	
	int is_dominanted;
	int ind_obj;

	initilize_dominance(dominance, solutions, size);
	for (int i =0; i< *size; i++){
		for (int j = 0; j< *size; j++){
			if (i != j){
				ind_obj = 0;
				is_dominanted = 1;
				//Check all objectives are less than other solution
				while ( (is_dominanted == 1) &&
						(ind_obj < dominance[j].sol->num_obj)){
					if (dominance[j].sol->obj_values[ind_obj] <= dominance[i].sol->obj_values[ind_obj]){
						is_dominanted = 1;
					}else{
						is_dominanted = 0;
					}
					ind_obj++;
				}
				if (is_dominanted == 1 ){
					is_dominanted = 0;
					// Check if an objective is less
					ind_obj = 0;
					while ( (is_dominanted == 0) &&
							(ind_obj < dominance[j].sol->num_obj)){						
						if (dominance[j].sol->obj_values[ind_obj] < dominance[i].sol->obj_values[ind_obj]){
							is_dominanted = 1;
						}
						ind_obj++;
					}
				}
				//i is dominated by j ????
				//Set value of how_many_is_dominated_by and set_dominated
				if (is_dominanted == 1){ //indicates that i is dominated by j
					// it shows how many solutions dominates i
					dominance[i].how_many_solutions_dominate_it = dominance[i].how_many_solutions_dominate_it + 1;
					// it indicates that i is dominated by j
					if  ( (dominance[j].max_dominated > *size) ||
							(dominance[j].max_dominated < 0) ){
						fatal_error("Index for ind_dominated is wrong. Check it! \n");
					}
					dominance[j].set_dominated[dominance[j].max_dominated] = i;
					dominance[j].max_dominated = dominance[j].max_dominated + 1;
				}
			}
		}
	}	
}

/** Initialize dominance struct
* dominance is a pointer that stores the application of dominance concept in solution
* solutions array of Solutions
* size number of solutions
*/
static void initilize_dominance(dominance_t *dominance, const solution_t *solutions, 
	const int *size){	
	for (int s = 0; s < *size; s++){
		dominance[s].sol = &solutions[s];
		dominance[s].max_dominated = 0;
		dominance[s].set_dominated = Malloc(int, *size);
		dominance[s].how_many_solutions_dominate_it = IS_DOMINANTED;
	}
}


