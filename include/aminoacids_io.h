#ifndef OLD_AMINOACIDS_IO_H
#define OLD_AMINOACIDS_IO_H

#include "aminoacids_types.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
primary_seq_t *_load_amino_seq(const char *file_name_protein);

#endif
