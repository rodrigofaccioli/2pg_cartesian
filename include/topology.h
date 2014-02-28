#ifndef OLD_TOPOLOGY_H
#define OLD_TOPOLOGY_H

#include "topology_types.h"
#include "protein_type.h"
#include "pdb_types.h"
#include "enums.h"

top_global_t *allocateTop_Global(const int *numres,
	const int *numatom);

void  desAllocateTop_Global(top_global_t *top_aux);
void build_topology_population(protein_t *pop, const int *pop_size);
void build_topology_individual(protein_t *prot);
int get_number_hydrogen_backbone(const protein_t *prot, const int *res);
boolean_t is_hydrogen_backbone_Nitrogen(const char *atomname);
boolean_t is_fixed_atom(const int *atmnumber, const int *fixed_atoms, const int *num_fixed);
int get_number_atoms_backbone(const protein_t *prot, const int *numres);
boolean_t is_backbone_atom(const char *atomname);

#endif
