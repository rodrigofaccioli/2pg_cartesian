#ifndef OLD_LOAD_PARAMETERS_H
#define OLD_LOAD_PARAMETERS_H

#include "parameters_type.h"

void load_parameters_from_file(input_parameters_t *in_param,
		const char *conf_file_name);
void set_parameters_extended_chain(input_parameters_extended_chain_t *param,
		char *argv[]);
void deAllocateload_parameters(input_parameters_t *param);
void deAllocateload_parameters_extended_chain(
		input_parameters_extended_chain_t *param);

#endif
