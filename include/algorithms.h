#ifndef OLD_ALGORITHMS_H
#define OLD_ALGORITHMS_H

#include "parameters_type.h"
#include "aminoacids_types.h"
#include "solution_types.h"

void initialize_algorithm_execution(const primary_seq_t *primary_sequence_aux,
		const input_parameters_t *in_para_aux);
void update_execution_algorithms(const solution_t *solutions, const int *tag);
int get_choose_residue(const int *num_res_prot);
#endif
