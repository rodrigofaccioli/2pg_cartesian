#ifndef OLD_DOMINANCE_H
#define OLD_DOMINANCE_H

#include "dominance_type.h"

dominance_t * allocate_dominance(const int *size);
void desallocate_dominance(dominance_t *dominance, const int *size);
void set_dominance(dominance_t *dominance, const solution_t *solutions, const int *size);
void show_dominance(const dominance_t *dominance, const int *size);
void save_dominance(const dominance_t *dominance, const int *size);

static void initilize_dominance(dominance_t *dominance, const solution_t *solutions, 
	const int *size);


#endif