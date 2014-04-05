#include <stdlib.h>
#include <string.h>

#include "load_parameters.h"
#include "random_algorithm.h"
#include "messages.h"

int main(int argc, char *argv[]){
	input_parameters_t in_param;
	display_msg("Reading the configure file \n");
	load_parameters_from_file(&in_param,argv[1]);

	//call random algorithm
	random_algorithm(&in_param);

	deAllocateload_parameters(&in_param);

	display_msg("Done Random Algorithm  !!! \n");
	return 0;
}