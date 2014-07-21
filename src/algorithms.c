#include <stdio.h>
#include <string.h>

#include "algorithms.h"
#include "algorithms_types.h"
#include "messages.h"
#include "defines.h"
#include "populationio.h"
#include "solutionio.h"
#include "objective.h"
#include "randomlib.h"
#include "rotation.h"
#include "math_owner.h"
#include "protein.h"
#include "topology.h"
#include "pdbatom.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

//#include "randomlib.h"

/* Starting declaration of global variable for algorithms */

/* These variables are used frequently. Instead of extern declaration, they were
 * declared as below because I would like to get more control for them.
 */
static const primary_seq_t *primary_sequence;
static const input_parameters_t *in_para;
static char pop_file_name[MAX_FILE_NAME];
static char fitness_file_name[MAX_FILE_NAME];
static char pop_file_name_GRO[MAX_FILE_NAME];
/* End of declaration of global variable for algorithms */

/** Initialization for executing of Algorithms */
_2PG_CARTESIAN_EXPORT
void initialize_algorithm_execution(const primary_seq_t *primary_sequence_aux,
		const input_parameters_t *in_para_aux){	
	primary_sequence = primary_sequence_aux;	
	in_para = in_para_aux;
}

static void set_objective_file_name(const int *fit, const int *generation){
    char fitness_name[MAX_RANDOM_STRING];
    char sger[MAX_RANDOM_STRING];
    type_fitness_energies2str(fitness_name, &in_para->fitness_energies[*fit]);
    sprintf(sger, "%d", *generation);
    strcat(fitness_name,"_");
    strcat(fitness_name,sger);
    sprintf(fitness_file_name,"%s.fit",fitness_name );
}

void build_fitness_files(const solution_t *solutions, const int *generation,
        const int *pop_size){
    int f;
    for (f =0; f < in_para->number_fitness; f++){
        set_objective_file_name(&f,generation);
        save_solution_file(in_para->path_local_execute, fitness_file_name, &f, 
            solutions, pop_size, generation, in_para);
    }
}

/** Building the name of population file  */
static void set_population_file_name(const int *tag){
    sprintf(pop_file_name,"pop_%d.pdb",*tag);
}


/** Updates the execution of Algorithms */
void update_execution_algorithms(const solution_t *solutions, const int *tag){
	const protein_t *population_aux;
    //Building population file name
    set_population_file_name(tag);

    //Saving population
    population_aux = (protein_t*) solutions[0].representation;
    save_population_file(population_aux, in_para->path_local_execute, pop_file_name, 
        &in_para->size_population);
    //Saving fitness values    
    build_fitness_files(solutions, tag, &in_para->size_population);
}

/**
*/
_2PG_CARTESIAN_EXPORT
int get_choose_residue(const int *num_res_prot){
    int max_residue_choose;
    int num_residue_choose;

    num_residue_choose = 0;    
    max_residue_choose = *num_res_prot + 1;
    while (num_residue_choose == 0){
        num_residue_choose = _get_int_random_number(&max_residue_choose);
    }
    return num_residue_choose;    
}

/** Implementation of 1 point crossover
* p_new represents the new solution
* p1 means the solution 1 that will build s_new
* p2 means the solution 2 that will build s_new
*/
static void crossover_1_part(protein_t *p_new, const protein_t *p1, 
    const protein_t *p2){
    int res_aux, res_ini;
    res_aux = 0;
    res_ini = 1;

    //getting a random residue
    res_aux = get_choose_residue(&p1->p_topol->numres);
    //coping from solution s1
    copy_protein_atoms_by_residues(p_new, &res_ini, &res_aux, p1);
    //coping from solution s2
    res_ini = res_aux + 1;
    if (res_ini <= p2->p_topol->numres){
        res_aux = p2->p_topol->numres;
        copy_protein_atoms_by_residues(p_new, &res_ini, &res_aux, p2);
    }
    build_topology_individual(p_new);
}

/** Applies the crossover parents.
* p_new is based on one of its parents
* When choose = 0 is considered father 1. Otherwise, father 2
*/
void crossover_parents(protein_t *p_new, const protein_t *p1,
        const protein_t *p2){
    int choose;
    int aux = 2;
    choose = _get_int_random_number(&aux);
    if (choose == 0){        
        copy_protein_atoms(p_new,p1);
    }else{
        copy_protein_atoms(p_new,p2);
    }
}

/** Applies the crossover operator
* ind_new is an individual of new population
* prot_1 is first indiviual
* prot_2 is second indiviual
*/
void apply_crossover(protein_t *ind_new, const protein_t *prot_1, 
    const protein_t * prot_2, type_crossoers_t *crossovers){
    if (crossovers[0] != crossoer_none){
        int choose;
        int aux = 2;
        choose = _get_int_random_number(&aux);
        if (choose == 0){
            crossover_parents(ind_new, prot_1, prot_2);
        }else{
            crossover_1_part(ind_new, prot_1, prot_2);
        }
    }else{
        crossover_parents(ind_new, prot_1, prot_2);
    }
}

/** Applies the mutation operator
* ind_new is an individual of new population
* in_para is the input parameter
*/
void apply_mutation(protein_t *ind_new, const input_parameters_t *in_para){
    float angle_radians, angle_degree, rate;
    int num_residue_choose;
    int what_rotation;
    int max_kind_of_rotation = 4;  

    rate = _get_float();
    num_residue_choose = 0;
    if (rate < in_para->individual_mutation_rate){
        for (int number_rotations = 1; number_rotations <= 
            in_para->how_many_rotations; number_rotations++){
            //Choose a residue
            num_residue_choose = get_choose_residue(&ind_new->p_topol->numres);
            //Obtaing kind of rotation
            what_rotation = _get_int_random_number(&max_kind_of_rotation);        
            //Appling the rotation 
            if (what_rotation == 0){
                //Obtaing a random degree angule 
                angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_psi, 
                &in_para->max_angle_mutation_psi);
                //Cast to radians
                angle_radians = degree2radians(&angle_degree);
                rotation_psi(ind_new, &num_residue_choose, &angle_radians);
            }else if (what_rotation == 1){
                //Obtaing a random degree angule 
                angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_phi, 
                &in_para->max_angle_mutation_phi);
                //Cast to radians
                angle_radians = degree2radians(&angle_degree);
                rotation_phi(ind_new, &num_residue_choose, &angle_radians);            
            }else if (what_rotation == 2){
                //Obtaing a random degree angule 
                angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_omega, 
                &in_para->max_angle_mutation_omega);
                //Cast to radians
                angle_radians = degree2radians(&angle_degree);
                rotation_omega(ind_new, &num_residue_choose, &angle_radians);            
            }else if (what_rotation == 3){
                int chi = 0;
                int max_chi = -1;
                char *res_name;
                res_name = Malloc(char, 4);
                get_res_name_from_res_num(res_name, &num_residue_choose, 
                    ind_new->p_atoms, &ind_new->p_topol->numatom);
                max_chi = get_number_chi(res_name);                        
                if (max_chi > 0){
                    //Choose a chi of residue. It must be started 1 untill number of residue
                    if (max_chi == 1){
                        chi = 1;
                    }else{
                        while (chi == 0){                
                            chi = _get_int_random_number(&max_chi);
                        }                    
                    }                
                    //Obtaing a random degree angule 
                    angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_side_chain, 
                    &in_para->max_angle_mutation_side_chain);
                    //Cast to radians
                    angle_radians = degree2radians(&angle_degree);
                    rotation_chi(ind_new, &num_residue_choose, &chi, &angle_radians);
                }
                free(res_name);
            }
        }
    }
}

/** Returns what generation want to start
*/
int get_started_generation(const int *start){
    return *start+1;
}
