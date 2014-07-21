#ifndef OLD_TOPOLOGY_H
#define OLD_TOPOLOGY_H

#include "topology_types.h"
#include "protein_type.h"
#include "pdb_types.h"
#include "enums.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
top_global_t *allocateTop_Global(const int *numres,
	const int *numatom);

_2PG_CARTESIAN_EXPORT
void  desAllocateTop_Global(top_global_t *top_aux);

_2PG_CARTESIAN_EXPORT
void build_topology_population(protein_t *pop, const int *pop_size);

_2PG_CARTESIAN_EXPORT
void build_topology_individual(protein_t *prot);

_2PG_CARTESIAN_EXPORT
int get_number_hydrogen_backbone(const protein_t *prot, const int *res);

_2PG_CARTESIAN_EXPORT
boolean_t is_hydrogen_backbone_Nitrogen(const char *atomname);

_2PG_CARTESIAN_EXPORT
boolean_t is_fixed_atom(const int *atmnumber, const int *fixed_atoms, const int *num_fixed);

_2PG_CARTESIAN_EXPORT
int get_number_atoms_backbone(const protein_t *prot, const int *numres);

_2PG_CARTESIAN_EXPORT
boolean_t is_backbone_atom(const char *atomname);

_2PG_CARTESIAN_EXPORT
int get_number_chi(const char *res_name);

_2PG_CARTESIAN_EXPORT
void rename_oxygen_c_terminal(pdb_atom_t *atoms,
		const int *res_num, const int *num_atom);

#endif
