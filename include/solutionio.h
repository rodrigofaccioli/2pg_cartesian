#ifndef OLD_SOLUTIONIO_H
#define OLD_SOLUTIONIO_H

#include "solution_types.h"
#include "parameters_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

void save_solution_file(const char *path, const char *file_name,
		const int *fit, const solution_t *solutions, const int *pop_size,
		const int *tag, const input_parameters_t *in_para);
_2PG_CARTESIAN_EXPORT
solution_t * loading_file_solutions(int *num_solutions_r, 
	int *numobj_r, 	const char *path_file_name);

#endif
