#ifndef OLD_TOPOLOGYIO_H
#define OLD_TOPOLOGYIO_H

#include "topology_types.h"
#include "protein_type.h"

void save_topology_protein(const top_global_t *top, const char *path, 
	const char *filename);
void save_topology_population(const protein_t *pop, const int *popsize, 
	const char *path, const char *prefix);

void _create_fasta_pdb(const char *prot_name, const char *chain_name,
		const char *prot_seq, const char *file_name_protein);

static boolean_t _check_pdb_fasta_file(char *line);
#endif
