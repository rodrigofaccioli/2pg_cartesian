#ifndef OLD_EA_MONO_H
#define OLD_EA_MONO_H

#include "parameters_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
int ea_mono(const input_parameters_t *in_para);

#endif
