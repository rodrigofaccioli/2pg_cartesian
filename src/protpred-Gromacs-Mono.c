#include <stdlib.h>
#include <string.h>

#include "load_parameters.h"
#include "ea_mono.h"
#include "messages.h"

int main(int argc, char *argv[]){
	input_parameters_t in_param;
	display_msg("Reading the configure file \n");
	load_parameters_from_file(&in_param,argv[1]);

	if (in_param.number_fitness == 1){
		ea_mono(&in_param);
	}else{
		fatal_error("In mono-objective, NumberObjective option must be 1. Please, check it. \n");
	}

	deAllocateload_parameters(&in_param);

	display_msg("Done EA Mono !!! \n");
	return 0;
}
