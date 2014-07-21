#ifndef OLD_PROTEIN_H
#define OLD_PROTEIN_H

#include "protein_type.h"
#include "solution_types.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
protein_t * allocateProtein(const int *size);

_2PG_CARTESIAN_EXPORT
void desallocateProtein(protein_t *pop, const int *inPopSize);

_2PG_CARTESIAN_EXPORT
void set_proteins2solutions(solution_t *sol, protein_t *pop, const int *pop_size);

_2PG_CARTESIAN_EXPORT
void copy_protein(protein_t *p_dest, const protein_t *p_source);

_2PG_CARTESIAN_EXPORT
void copy_protein_population(protein_t *pop_dest, const protein_t *pop_source, 
	const int *popsize);

_2PG_CARTESIAN_EXPORT
void copy_protein_atoms(protein_t *p_dest, const protein_t *p_source);

_2PG_CARTESIAN_EXPORT
void copy_protein_atoms_by_residues(protein_t *p_dest, const int *res_num_ini, 
	const int *res_num_end, const protein_t *p_source);

_2PG_CARTESIAN_EXPORT
void copy_protein_population_atoms(protein_t *pop_dest, const protein_t *pop_source, 
	const int *popsize);

_2PG_CARTESIAN_EXPORT
void initialize_protein_atoms(protein_t *protein);

_2PG_CARTESIAN_EXPORT
void initialize_protein_population_atoms(protein_t *protein, const int *popsize);
#endif
