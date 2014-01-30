#include "protein.h"

void write_initial_population_file(const char *path, const char *file_name,
		protein ** pop, const int *pop_size);
void _save_population_file(const char *path, const char *file_name,
		protein ** pop, const int *pop_size);
void _save_protein_path_file(const char *path_file_name, const protein *prot);

static void start_ind(FILE *pop_file);
static void finish_ind(FILE *pop_file);
static void write_individual(FILE *pop_file,protein **pop,
		const int *ind);
static void write_protein(FILE *pop_file, const protein *prot);
