#ifndef OLD_AMINOACIDS_IO_H
#define OLD_AMINOACIDS_IO_H

#include "aminoacids_types.h"

#ifdef __cplusplus
extern "C"
{
#endif

primary_seq_t *_load_amino_seq(const char *file_name_protein);

#ifdef __cplusplus
}
#endif

#endif