#include "topology_types.h"
#include "vector_types.h"
#include "pdb_types.h"


float _compute_diehdral_angle(const own_vector_t *a1,
		const own_vector_t *a2,const own_vector_t *a3,	const own_vector_t *a4);

float _compute_phi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global);

float _compute_psi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global);

float  _compute_side_chains_angles(const int *r, const int *chi,
		const pdb_atom_t *pdb_atoms, const top_global_t *top_global);

void _type_of_diedhral_angle2str(char *str,
		const type_dihedral_angles_t *type_dihedral);

static int _get_atom_index_from_top_global_dihedral_side_chain_t(
		int atm_opt, const int *r, 	const int *chi,
		const top_global_t *top_global);