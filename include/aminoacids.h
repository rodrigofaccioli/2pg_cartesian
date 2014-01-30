#ifndef OLD_AMINOACIDS_H
#define OLD_AMINOACIDS_H

#include "aminoacids_types.h"

primary_seq_t* allocate_primary_seq(const int *num_res);
void desallocate_primary_seq(primary_seq_t* seq);

#endif