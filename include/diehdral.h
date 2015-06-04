#ifndef OLD_DIEHDRAL_H
#define OLD_DIEHDRAL_H

#include "topology_types.h"
#include "vector_types.h"
#include "pdb_types.h"

float compute_diehdral_angle(const own_vector_t *a1,
		const own_vector_t *a2,const own_vector_t *a3,	
		const own_vector_t *a4);
float compute_phi_residue(pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top);
float compute_psi_residue(pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top);
float compute_omega_residue(pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top);
void compute_chi_residue(float *chi, int *num_chi, pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top);

#endif