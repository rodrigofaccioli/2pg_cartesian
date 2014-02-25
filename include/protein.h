#ifndef OLD_PROTEIN_H
#define OLD_PROTEIN_H

#include "protein_type.h"
#include "solution_types.h"

protein_t * allocateProtein(const int *size);
void desallocateProtein(protein_t *pop, const int *inPopSize);
void set_proteins2solutions(solution_t *sol, protein_t *pop, const int *pop_size);

#endif
