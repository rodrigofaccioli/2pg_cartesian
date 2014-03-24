#ifndef OLD_TOPOLOGY_LIB_H
#define OLD_TOPOLOGY_LIB_H


#include "topology_types.h"
#include "pdb_types.h"
#include "protein_type.h"

int get_atom_index_by_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname, const int *num_atom);
boolean_t atom_name_exists_in_resnum(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname, const int *num_atom);
void set_fixed_moved_atoms_side_chains_SER(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_CYS(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_THR(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_VAL(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_ASP(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_ASN(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_LEU(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_PRO(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_PHE(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_HIS(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_TYR(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_TRP(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_MET(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_GLN(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_GLU(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_LYS(protein_t *prot, 
	const int *res_num, const int *chi);
void set_fixed_moved_atoms_side_chains_ARG(protein_t *prot, 
	const int *res_num, const int *chi);
#endif


