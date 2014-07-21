#ifndef OLD_ROTATION_TYPE_H
#define OLD_ROTATION_TYPE_H

#include "protein_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
void rotation_psi(protein_t *prot, const int *num_res_first, const float *angle);
_2PG_CARTESIAN_EXPORT
void rotation_phi(protein_t *prot, const int *num_res_first, const float *angle);
_2PG_CARTESIAN_EXPORT
void rotation_omega(protein_t *prot, const int *num_res_first, const float *angle);
_2PG_CARTESIAN_EXPORT
void rotation_chi(protein_t *prot, const int *num_res_first, const int *chi,
	const float *angle);

#endif
