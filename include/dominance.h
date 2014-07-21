#ifndef OLD_DOMINANCE_H
#define OLD_DOMINANCE_H

#include "dominance_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
dominance_t * allocate_dominance(const int *size);
_2PG_CARTESIAN_EXPORT
void desallocate_dominance(dominance_t *dominance, const int *size);
_2PG_CARTESIAN_EXPORT
void set_dominance(dominance_t *dominance, const solution_t *solutions, const int *size);
_2PG_CARTESIAN_EXPORT
void show_dominance(const dominance_t *dominance, const int *size);
_2PG_CARTESIAN_EXPORT
void save_dominance(const dominance_t *dominance, const int *size);

static void initilize_dominance(dominance_t *dominance, const solution_t *solutions, 
	const int *size);


#endif
