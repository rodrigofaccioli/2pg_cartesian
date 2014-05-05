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
	if (p_dest->p_atoms == NULL){
		p_dest->p_atoms = allocate_pdbatom(&p_source->p_topol->numatom);
	}
	for (int a = 0; a < p_source->p_topol->numatom; a++){
		copy_pdb_atom(&p_dest->p_atoms[a], &p_source->p_atoms[a]);
	}
	if (p_dest->p_topol == NULL){
		p_dest->p_topol = allocateTop_Global(&p_source->p_topol->numres, 
			&p_source->p_topol->numatom);
	}
	build_topology_individual(p_dest);
}

/** Copies p_source to p_dest at reproduce*/
void copy_protein_at_reproduce(protein_t *p_dest, const protein_t *p_source){
	/* Checking if is necessary desalocate atoms and topology because topology, 
	* and consequently number of atoms, can be different between p_dest and p_souce.
	* Important: The difference betwwen dest and source occur in reprodution since
	* conformations od protein can be getted by different methods and it can change number
	* of atoms. Example, a method should add protonation state while others not.
	*/
	if (p_dest->p_atoms != NULL){
		desAllocate_pdbatom(p_dest->p_atoms);
		p_dest->p_atoms = NULL;
	}
	if (p_dest->p_topol != NULL){
		desAllocateTop_Global(p_dest->p_topol);
		p_dest->p_topol = NULL;
	}
	copy_protein(p_dest, p_source);
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

/** Copies atoms of p_source to p_dest based on residue
* p_dest means protein that will receive the atoms
* res_num_ini is number of residue starting the copy of atoms
* res_num_end is number of residue finishing the copy of atoms
* p_source means protein where the atoms are
*
* Important: res_ini and res_end are considered the number of
* residue. Therefore, it is not index of residue. 
*/
void copy_protein_atoms_by_residues(protein_t *p_dest, const int *res_num_ini, 
	const int *res_num_end, const protein_t *p_source){
	if (*res_num_ini <= 0){
		fatal_error("In copy_protein_atoms_by_residues function the value of res_num_ini must be more than zero!");
	}
	for (int r = *res_num_ini; r <= *res_num_end; r++ ){
		//Coping atoms
		for (int a = p_source->p_topol->range_atoms[r-1].first_atom; a < p_source->p_topol->range_atoms[r].last_atom; a++){			
			copy_pdb_atom(&p_dest->p_atoms[a-1], &p_source->p_atoms[a-1]);
		}
	}
}

/** Initializes the atoms of protein */
void initialize_protein_atoms(protein_t *protein){
	//set 0 to coordenates of protein
	for (int a = 0; a < protein->p_topol->numatom; a++){			
		protein->p_atoms[a].coord.x = 0.00;
		protein->p_atoms[a].coord.y = 0.00;
		protein->p_atoms[a].coord.z = 0.00;
	}
}

/** Initializes the atoms of each individual of protein */
void initialize_protein_population_atoms(protein_t *protein, const int *popsize){
	for (int i = 0; i < *popsize; i++){
		initialize_protein_atoms(protein);
	}
}
