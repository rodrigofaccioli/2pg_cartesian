#ifndef OLD_SOLUTIONIO_H
#define OLD_SOLUTIONIO_H

#include "solution_types.h"
#include "parameters_type.h"

void save_solution_file(const char *path, const char *file_name,
		const int *fit, const solution_t *solutions, const int *pop_size,
		const int *tag, const input_parameters_t *in_para);
solution_t * loading_file_solutions(int *num_solutions_r, 
	int *numobj_r, 	const char *path_file_name);
solution_t * loading_file_solutions_kind_objectives(int *num_solutions_r, 
	int *numobj_r, 	const char *path_file_name);

#endif
