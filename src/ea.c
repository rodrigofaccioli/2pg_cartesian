#include <stdio.h>
#include <string.h>

#include "ea.h"
#include "ea_types.h"
#include "randomlib.h"
#include "messages.h"
#include "defines.h"
#include "pdbio.h"
#include "topology.h"
#include "topologylib.h"

/* Starting declaration of global variable for ea */

/* These variables are used frequently. Instead of extern declaration, they were
 * declared as below because I would like to get more control for them.
 */
static const primary_seq_t *primary_sequence = NULL;
static const input_parameters_t *in_para;
static char pop_file_name[MAX_FILE_NAME];
static char fitness_file_name[MAX_FILE_NAME];
static char pop_file_name_GRO[MAX_FILE_NAME];
/* End of declaration of global variable for ea */

void ea_initialize(const primary_seq_t *primary_sequence_aux,
		const input_parameters_t *in_para_aux){
	/* This function prepares to execution of Evolutionary Algorithms.
	 * Firstly, it is loading database which is used to execute the mutation
	 * genetic operators.
	 */
	display_msg("Initializing global variables \n");
	primary_sequence = primary_sequence_aux;	
	in_para = in_para_aux;
}

