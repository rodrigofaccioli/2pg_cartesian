#ifndef OLD_TOPOLOGY_H
#define OLD_TOPOLOGY_H

#include "protein.h"
#include "pdbatom.h"
#include "topology_types.h"
#include "parameters_type.h"

top_global_t *allocateTop_Global(const int *numatom, const int *numres,
		const int *num_bond_angles, const int *num_side_chains,
		const int *number_dihedrals_type);

void  desAllocateTop_Global(top_global_t *top_global);

#endif
