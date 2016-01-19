#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "load_parameters.h"
#include "LoadConfig.h"
#include "messages.h"
#include "string_owner.h"
#include "math_owner.h"
#include "mutations.h"

type_mutations_t str2type_mutation(char *name_mutation)
{
	if (strcmp(name_mutation,"general_rotation") == 0)
    {
	    return general_rotation;
	}
    else if (strcmp(name_mutation,"dm_refinement") == 0)
    {
		return dm_refinement;
	}
    else if (strcmp(name_mutation,"none") == 0)
    {
        return none;
    }
    else
    {
		char msg[1024];
		sprintf(msg,"The option %s did not find at type_mutations_t. Please, check it!", name_mutation);
		fatal_error(msg);
	}
}

static void set_parameter_mutations(input_parameters_t *param, char *mutations_parameters)
{
    char *tokens;
    int count=1; /*initialized with the value "1" , because we can assume that at least the "None" value in the "Mutation_Operator" - check the file "configuracao.conf" */

    while((tokens = strpbrk(mutations_parameters, ",")) != NULL)
    {
        count++;
        tokens++;
    }

    in_para->mutations = Malloc(type_mutations_t, count);
    int i=0;
    tokens = strtok(mutations_parameters, " ,");

    while(tokens != NULL)
    {
         in_para->mutations[i] = str2type_mutation(tokens);
         i++;

         if(count>1 && (strcmp(tokens, "None") == 0))
         {
            char msg[300];
            sprintf(msg, "The number of mutations are %d. \"None\" along with other mutations in parameters. Check it!", count);
            fatal_error(msg);

         }

         tokens = strtok(NULL, " ,");
    }
}

void general_rotation(protein_t *ind_new, const input_parameters_t *in_para)
{
    float angle_radians, rate;
    int num_residue_choose;
    int what_rotation;
    int max_kind_of_rotation = 4;

    rate = _get_float();
    num_residue_choose = 0;

    if (rate < in_para->individual_mutation_rate)
    {
        for (int number_rotations = 1; number_rotations <= in_para->how_many_rotations; number_rotations++)
        {
            //Choose a residue
            num_residue_choose = get_choose_residue(&ind_new->p_topol->numres);
            //Obtaing kind of rotation
            what_rotation = _get_int_random_number(&max_kind_of_rotation);
            //Appling the rotation
            switch(what_rotation)
            {
                case 0:
                        angle_radians = _get_float_random_interval(&in_para->min_angle_mutation_psi,
                        &in_para->max_angle_mutation_psi);
                        rotation_psi(ind_new, &num_residue_choose, &angle_radians);
                        break;
                case 1:
                        angle_radians = _get_float_random_interval(&in_para->min_angle_mutation_phi,
                        &in_para->max_angle_mutation_phi);
                        rotation_phi(ind_new, &num_residue_choose, &angle_radians);
                        break;
                case 2:
                        //Obtaing a random degree angule
                        angle_radians = _get_float_random_interval(&in_para->min_angle_mutation_omega,
                        &in_para->max_angle_mutation_omega);
                        rotation_omega(ind_new, &num_residue_choose, &angle_radians);
                        break;
                case 3:
                        int chi = 0;
                        int max_chi = -1;
                        char *res_name;
                        res_name = Malloc(char, 4);
                        get_res_name_from_res_num(res_name, &num_residue_choose,
                            ind_new->p_atoms, &ind_new->p_topol->numatom);
                        max_chi = get_number_chi(res_name);

                        switch(max_chi)
                        {
                            case 1:
                                    chi = 1;
                                    break;
                            case 0:
                                    chi = _get_int_random_number(&max_chi);
                                    break;
                            default:
                                    break;
                        }
                        free(res_name);
                        break;
                default:
                        break;
            }
        }
    }
}
