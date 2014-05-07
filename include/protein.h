#ifndef OLD_PROTEIN_H
#define OLD_PROTEIN_H

#include "protein_type.h"
#include "solution_types.h"

protein_t * allocateProtein(const int *size);
void desallocateProtein(protein_t *pop, const int *inPopSize);
void set_proteins2solutions(solution_t *sol, protein_t *pop, const int *pop_size);
void copy_protein(protein_t *p_dest, const protein_t *p_source);
void copy_protein_population(protein_t *pop_dest, const protein_t *pop_source, 
	const int *popsize);
void copy_protein_atoms(protein_t *p_dest, const protein_t *p_source);
void copy_protein_atoms_by_residues(protein_t *p_dest, const int *res_num_ini, 
	const int *res_num_end, const protein_t *p_source);
void copy_protein_population_atoms(protein_t *pop_dest, const protein_t *pop_source, 
	const int *popsize);
void initialize_protein_atoms(protein_t *protein);
void initialize_protein_population_atoms(protein_t *protein, const int *popsize);
#endif
