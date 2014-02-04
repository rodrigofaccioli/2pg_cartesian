#ifndef OLD_TOPOLOGY_H
#define OLD_TOPOLOGY_H

#include "topology_types.h"
#include "aminoacids_types.h"

top_global_t *allocateTop_Global(const primary_seq_t *primary_sequence,
	const int *numatom);

void  desAllocateTop_Global(top_global_t *top_global);

#endif
