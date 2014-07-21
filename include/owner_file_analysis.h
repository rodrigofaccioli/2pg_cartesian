#ifndef OLD_OWNER_FILE_H
#define OLD_OWNER_FILE_H

#include "analysis_types.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
owner_file_t* allocate_file_t(const int *num_files, const int *num_obj);
_2PG_CARTESIAN_EXPORT
void desalocate_file_t(owner_file_t*file_names, const int *num_files);
void copy_file_owner(owner_file_t* dest, const owner_file_t*source, const int *num);
const char *get_filename_ext(const char *filename);
int how_many_files_directory_by_extension(const char *path, const char *ext);
void insert_files_directory_by_extension(owner_file_t *file_names, const char *path, const char *ext);
double get_specific_objective_from_owner_file_t(const owner_file_t *sol, const int *obj);
double get_oposite_specific_objective_from_owner_file_t(const owner_file_t *sol, const int *obj);
double get_displayed_value_of_objective_from_owner_file_t(const owner_file_t *sol, 
	const int *index, const int *obj, 
	const type_fitness_energies_t *fitness_energies);
int compute_how_many_front_file_t(const owner_file_t * solutions, 
            const int *size, const int *front_ref);	
int compare_first_value(const void *x, const void *y);
void save_analysis_files(const owner_file_t *solutions_f, const int *size, const int *numobj, const type_fitness_energies_t *fitness_energies);
void save_analysis_files_no_objectives(const owner_file_t *solutions_f, const int *size, const int *numobj, const char *column_1, const char *column_2 );
int compare_front(const void *x, const void *y);
int compare_dominated(const void *x, const void *y);
void set_global_ranking(owner_file_t *file_names, const int *size);
void sorting_solutions_by_front_dominance(owner_file_t *file_names, const int *size, const int *numobj);
owner_file_t *loading_owner_file_solution(const int *num_solutions_r, const int *numobj_r, const char *path_file_name);

#endif
