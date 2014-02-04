#ifndef OLD_POPULATION_IO_H
#define OLD_POPULATION_IO_H

#include "protein_type.h"
#include "aminoacids_types.h"

void load_initial_population_file(protein_t *pop, const int *pop_size, const char *path, 
	const char *file_name, const primary_seq_t *primary_sequence);
void save_population_file(const protein_t *pop, const char *path, const char *file_name, 
	const int *num_model );
#endif