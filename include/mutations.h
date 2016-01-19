#ifndef OLD_MUTATIONS_H
#define OLD_MUTATIONS_H

#include "parameters_type.h"
#include "enums.h"


#ifdef __cplusplus
extern "C"
{
#endif

type_mutations_t str2type_mutation(char *name_mutation);
static void set_parameter_mutations(input_parameters_t *param, char *mutations_parameters);
void general_rotation(protein_t *ind_new, const input_parameters_t *in_para);

#ifdef __cplusplus
}
#endif

#endif
