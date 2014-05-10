#include <stdlib.h>
#include <string.h>

#include "dominance.h"
#include "load_parameters.h"
#include "messages.h"
#include "defines.h"
#include "solution.h"
#include "futil.h"
#include "string_owner.h"
#include "solutionio.h"
#include "ea_nsga2.h"

static dominance_t *dominance;
static solution_t  *solutions;
static ea_nsga2_t *nsga2_solutions_p; 
static int num_solutions;
static int num_obj;

void save_front(const ea_nsga2_t *nsga2_solutions, const int *size){
	FILE *d_file=NULL;
	char *file_name;
	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "front_output.txt");	
    
	d_file = open_file(file_name, fWRITE);

	fprintf(d_file, "#Identification\tFront\n");
	for (int s = 0; s < *size; s++){
		fprintf(d_file, "%i\t%i\n", nsga2_solutions[s].sol->ID, nsga2_solutions[s].front);
	}
	fclose(d_file);
	free(file_name);	
}

static int compare_front(const void *x, const void *y){
    int fx, fy;
    fx = ((ea_nsga2_t *)x)->front;
    fy = ((ea_nsga2_t *)y)->front;
    if (fx > fy){
        return 1;
    }else {
        return 0;
    }
}

int main(int argc, char *argv[]){
	input_parameters_t in_param;
	//display_msg("Reading the configure file \n");
	//load_parameters_from_file(&in_param,argv[1]);


	char *path_file_name;
	path_file_name = Malloc(char, MAX_PATH_FILE_NAME);
	strcpy(path_file_name, argv[1]);	

	//Getting Solutions
	solutions = loading_file_solutions(&num_solutions, &num_obj, path_file_name);
	in_param.size_population = num_solutions;
	in_param.number_fitness = num_obj;
	nsga2_solutions_p = allocate_nsga2(&in_param);	
	//Setting identification
	for (int i = 0; i < num_solutions; i++){
		nsga2_solutions_p[i].sol->ID = i+1;
	}

    //Setting dominance
	dominance = allocate_dominance(&num_solutions);
	set_dominance(dominance, solutions, &num_solutions);

    //Coping values of objective    
    set_nsga2_solution_in_solution(solutions, nsga2_solutions_p, &num_solutions);
	
    //Setting front based on dominance concept
    compute_fronts(nsga2_solutions_p, dominance, &num_solutions);

    //Sorting by front
    qsort(nsga2_solutions_p, num_solutions,  sizeof (ea_nsga2_t), compare_front);

    //Saving file
	save_front(nsga2_solutions_p, &num_solutions);

	desallocate_dominance(dominance, &num_solutions);
	desallocate_solution(solutions, &num_solutions);
	desallocate_solution_nsga2(nsga2_solutions_p, &num_solutions);
	free(path_file_name);

	display_msg("Done Computing Front !!! \n");

	return 0;
}