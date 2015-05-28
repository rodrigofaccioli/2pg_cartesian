#ifndef OLD_ROTATION_TYPE_H
#define OLD_ROTATION_TYPE_H

#include "protein_type.h"

void rotation_psi(protein_t *prot, const int *num_res_first, const float *angle);
void rotation_psi_residue(protein_t *prot, const int *num_res, const float *angle);
void rotation_phi(protein_t *prot, const int *num_res_first, const float *angle);
void rotation_phi_residue(protein_t *prot, const int *num_res, const float *angle);
void rotation_omega(protein_t *prot, const int *num_res_first, const float *angle);
void rotation_omega_residue(protein_t *prot, const int *num_res, const float *angle);
void rotation_chi(protein_t *prot, const int *num_res_first, const int *chi,
	const float *angle);

#endif