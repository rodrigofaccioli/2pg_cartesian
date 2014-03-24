#include <stdlib.h>
#include <string.h>

#include "load_parameters.h"
#include "mc_metropolis.h"
#include "messages.h"

int main(int argc, char *argv[]){
	input_parameters_t in_param;
	display_msg("Reading the configure file \n");
	load_parameters_from_file(&in_param,argv[1]);

	if (in_param.number_fitness == 1){
		mc_metropolis(&in_param);
	}else{
		fatal_error("In Monte Carlo Metropolis, NumberObjective option must be 1. Please, check it. \n");
	}

	deAllocateload_parameters(&in_param);

	display_msg("Done MC Metropolis !!! \n");
	return 0;
}