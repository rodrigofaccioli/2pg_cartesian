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

static dominance_t *dominance;
static solution_t  *solutions;
static int num_solutions;
static int num_obj;


int main(int argc, char *argv[]){
	input_parameters_t in_param;

	char *path_file_name;
	path_file_name = Malloc(char, MAX_PATH_FILE_NAME);
	strcpy(path_file_name, argv[1]);	

	solutions = loading_file_solutions(&num_solutions, &num_obj, path_file_name);

	dominance = allocate_dominance(&num_solutions);
	set_dominance(dominance, solutions, &num_solutions);
	show_dominance(dominance, &num_solutions);
	save_dominance(dominance, &num_solutions);
	
	desallocate_dominance(dominance, &num_solutions);
	desallocate_solution(solutions, &num_solutions);
	free(path_file_name);

	display_msg("Done Dominance !!! \n");

	return 0;
}
