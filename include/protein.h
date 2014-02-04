#ifndef OLD_PROTEIN_H
#define OLD_PROTEIN_H

#include "protein_type.h"

protein_t * allocateProtein(const int *size);
void desallocateProtein(protein_t *pop, const int *inPopSize);

#endif
