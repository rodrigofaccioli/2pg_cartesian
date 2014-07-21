#ifndef OLD_LOAD_PARAMETERS_H
#define OLD_LOAD_PARAMETERS_H

#include "parameters_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
void load_parameters_from_file(input_parameters_t *in_param,
		const char *conf_file_name);

_2PG_CARTESIAN_EXPORT
void deAllocateload_parameters(input_parameters_t *param);

#endif
