#ifndef OLD_OBJECTIVE_IO_H
#define OLD_OBJECTIVE_IO_H

#include "parameters_type.h"

void _save_objective_file(const char *path, const char *file_name,
		const int *fit, protein ** pop, const int *pop_size,
		const int *generation, const input_parameters_t *in_para);
#endif
