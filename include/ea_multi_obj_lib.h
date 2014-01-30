#ifndef OLD_EA_MULTI_OBJ_LIB_H
#define OLD_EA_MULTI_OBJ_LIB_H

#include "protein.h"
#include "z_matrix.h"
#include "parameters_type.h"
#include "topology.h"


void _fitness_gromacs(protein **population_p,
		const input_parameters_t *in_para, const top_global_t *top);
void build_random_individuals_multi(protein **p_new, const int *sizepop,
		const int *number_ind_select_reproduce);

#endif
