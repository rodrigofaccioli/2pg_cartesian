#ifndef OLD_PDBATOM_H
#define OLD_PDBATOM_H

#include "enums.h"
#include "pdb_types.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
pdb_atom_t** allocate_Population_pdb(const int *inPopSize, const int *numatom_by_model);
_2PG_CARTESIAN_EXPORT
void desAllocate_Population_pdb(pdb_atom_t** pdbatoms, const int *inPopSize);
_2PG_CARTESIAN_EXPORT
pdb_atom_t * allocate_pdbatom(const int *numatom);
_2PG_CARTESIAN_EXPORT
void desAllocate_pdbatom(pdb_atom_t *pdbatoms);
_2PG_CARTESIAN_EXPORT
void copy_pdb_atom(pdb_atom_t *dest, const pdb_atom_t *source);
_2PG_CARTESIAN_EXPORT
void set_pdb_atom_coordinates(pdb_atom_t *pdbatom,	char *atmname, char *resname,
		const char *chain_name, const int *resnum, 	const float *x,
		const float *y, const float *z, const int *index);
_2PG_CARTESIAN_EXPORT
void set_pdb_atom_generic_information(pdb_atom_t *pdbatom,
		char *atmname, char *resname,
		const char *chain_name, const int *resnum,  const int *index);
_2PG_CARTESIAN_EXPORT
void get_atom_name_from_atomid(char *atomname, const type_atoms_t *atomid);
_2PG_CARTESIAN_EXPORT
type_atoms_t get_atomid_from_atom_name(const  char *__atmname);
_2PG_CARTESIAN_EXPORT
const pdb_atom_t * search_pdb_atom_from_resnum_atomid(const pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom);
const pdb_atom_t * search_pdb_atom_from_resnum_atom_name(const pdb_atom_t *atoms,
_2PG_CARTESIAN_EXPORT
		const int *res_num, const char *atomname,	const int *num_atom);
_2PG_CARTESIAN_EXPORT
const pdb_atom_t * get_pdb_atom_from_resnum_atomid(const pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom);
_2PG_CARTESIAN_EXPORT
const pdb_atom_t * get_pdb_atom_from_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname,	const int *num_atom);
_2PG_CARTESIAN_EXPORT
int get_last_number_res_from_atom(const pdb_atom_t *atoms, const int *num_atom);
_2PG_CARTESIAN_EXPORT
int get_number_residues_from_atom(const pdb_atom_t *atoms, const int *num_atom);
_2PG_CARTESIAN_EXPORT
void get_res_name_from_res_num(char *res_name, const int *num_res,
		const pdb_atom_t *atoms, const int *num_atom);
_2PG_CARTESIAN_EXPORT
void renumerate_residue_number(pdb_atom_t *atoms, const int *num_atom);
_2PG_CARTESIAN_EXPORT
void rename_atom(pdb_atom_t *atoms, const char *name, const char *name_new,
		const int *res_num, const type_atoms_t *atomid,
		const type_atoms_t *atomid_new, const int *num_atom);
static boolean_t is_residue_number_ok(pdb_atom_t *atoms);
static pdb_atom_t * search_pdb_atom_from_resnum_atomid_alow_change(pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom);
#endif
