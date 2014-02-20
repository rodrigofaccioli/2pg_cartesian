#include <stdio.h>
#include <string.h>

#include "protein.h"
#include "defines.h"
#include "messages.h"
#include "pdbatom.h"
#include "topology.h"


protein_t * allocateProtein(const int *size){
	/** Allocates a population of protein_t as array.
	 */
	protein_t *protein_aux;
	protein_aux = Malloc(protein_t,*size);
	for (int i = 0; i < *size; i++){
		protein_aux[i].p_atoms = NULL;
		protein_aux[i].p_topol = NULL;
	}
    return protein_aux;
}

void desallocateProtein(protein_t *pop, const int *inPopSize){
	for (int i =0; i < *inPopSize; i++){
		if (pop[i].p_atoms != NULL){
			desAllocate_pdbatom(pop[i].p_atoms);
		}
		if (pop[i].p_topol != NULL){
			desAllocateTop_Global(pop[i].p_topol);			
		}
	}
	free(pop);
}
