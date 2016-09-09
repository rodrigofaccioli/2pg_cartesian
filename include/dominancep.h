#ifndef OLD_DOMINANCEP_H
#define OLD_DOMINANCEP_H

#include "dominancep_type.h"

dominancep_t * allocate_dominancep(const int *size);
void desallocate_dominancep(dominancep_t *dominancep, const int *size);
void set_dominancep(dominancep_t *dominancep, const solution_t *solutions, const int *size);
void show_dominancep(const dominancep_t *dominancep, const int *size);
void save_dominancep(const dominancep_t *dominancep, const int *size);

static void initialize_dominancep(dominancep_t *dominancep, const solution_t *solutions, 
	const int *size);


#endif
