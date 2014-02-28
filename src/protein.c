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
	/*
	for (int i =0; i < *inPopSize; i++){
		if (pop[i].p_atoms != NULL){
			desAllocate_pdbatom(pop[i].p_atoms);
		}
		if (pop[i].p_topol != NULL){
			desAllocateTop_Global(pop[i].p_topol);			
		}
	}*/
	free(pop);
}

/** Set a population of protein_t to solution
* sol means the solutions that will be reference to individual of population
* pop is the population which has the protein conformations
* pop_size size of population
*/
void set_proteins2solutions(solution_t *sol, protein_t *pop, const int *pop_size){
	for (int i = 0; i < *pop_size; i++){
		sol[i].representation = &pop[i];
	}
}

/** Copies p_source to p_dest */
void copy_protein(protein_t *p_dest, const protein_t *p_source){
	//Coping atoms
	for (int a = 0; a < p_source->p_topol->numatom; a++){			
		copy_pdb_atom(&p_dest->p_atoms[a], &p_source->p_atoms[a]);
	}		
	//Building Topology of p_dest
	build_topology_individual(p_dest);
}

/** Copies each individual of pop_source to pop_dest */
void copy_protein_population(protein_t *pop_dest, const protein_t *pop_source, 
	const int *popsize){
	for (int i = 0; i < *popsize; i++){
		copy_protein(&pop_dest[i], &pop_source[i]);
	}
}

/** Copies atoms of p_source to p_dest */
void copy_protein_atoms(protein_t *p_dest, const protein_t *p_source){
	//Coping atoms
	for (int a = 0; a < p_source->p_topol->numatom; a++){			
		copy_pdb_atom(&p_dest->p_atoms[a], &p_source->p_atoms[a]);
	}
}

/** Copies atoms of each individual of pop_source to pop_dest */
void copy_protein_population_atoms(protein_t *pop_dest, const protein_t *pop_source, 
	const int *popsize){
	for (int i = 0; i < *popsize; i++){
		copy_protein_atoms(&pop_dest[i], &pop_source[i]);
	}
}


