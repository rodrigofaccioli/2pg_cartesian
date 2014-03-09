#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "enums.h"
#include "ea_mono.h"
#include "protein.h"
#include "topology.h"
#include "pdbio.h"
#include "pdbatom.h"
#include "messages.h"
#include "algorithms.h"
#include "string_owner.h"
#include "futil.h"
#include "random_number_gsl.h"
#include "aminoacids.h"
#include "aminoacids_io.h"
#include "populationio.h"
#include "topology.h"
#include "topologyio.h"
#include "rotation.h"
#include "math_owner.h"
#include "solution.h"
#include "gromacs_objectives.h"
#include "algorithms.h"
#include "randomlib.h"

 /** Implementation of tournament for chosing individual to reproduce 
 * solutions represents the population of solutions
 * in_para represents the parameters
 * output is the reference of solution has earned the tournament
 */
static int tournament(const solution_t *solutions, const int *popsize){
    int ind_1, ind_2;
    int max_random;
    max_random = *popsize -1;
    ind_1 = _get_int_random_number(&max_random);
    ind_2 = _get_int_random_number(&max_random);
    if (solutions[ind_1].obj_values[0] < solutions[ind_2].obj_values[0]){
        return ind_1;
    }else{
        return ind_2;
    }
}

/** Selects the individuals which will be reproduced 
* index_solutions_to_reproduce index of individuals were selected to reproduce
* num_ind_reproduce number the individuals for reprodution
* solutions population of solutions
* popsize size of population
*/
static void select_individual_mono(int *index_solutions_to_reproduce, 
        const int *num_ind_reproduce, const solution_t *solutions, const int *popsize){    
    for (int i = 0; i < *num_ind_reproduce;i++){
       index_solutions_to_reproduce[i] = tournament(solutions,popsize);
    }

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
    const protein_t * prot_2){

    crossover_parents(ind_new, prot_1, prot_2);
}

/** Applies the mutation operator
* ind_new is an individual of new population
* in_para is the input parameter
*/
void apply_mutation(protein_t *ind_new, const input_parameters_t *in_para){
    float angle_radians, angle_degree, rate;
    int num_residue_choose;
    int what_rotation;
    int max_kind_of_rotation = 3;  

    rate = _get_float();
    num_residue_choose = 0;
    if (rate < in_para->individual_mutation_rate){
        //Choose a residue. It must be started 1 untill number of residue
        while (num_residue_choose == 0){
            num_residue_choose = _get_int_random_number(&ind_new->p_topol->numres);
        }
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
        }
    }

}

/** Applies the reprodution in proteins that are into solutions
* pop_new indicates the new population
* solutions contain current protein through the field called representation 
* index_solutions_to_reproduce shows the index of selected solutions
* in_para stores the input parameters
* size_to_reproduce indicates the size of population to reproduce. It can be
                    equals with population size
*/
void reproduce_protein(protein_t *pop_new, const solution_t *solutions, 
    const int *index_solutions_to_reproduce, 
    const input_parameters_t *in_para, const int *size_to_reproduce){

    protein_t *prot_aux_1, *prot_aux_2;
    int index_aux, ind_1, ind_2, max_random;

    max_random = in_para->size_population -1;
    for (int i = 0; i < *size_to_reproduce; i++){
        //Obtaining two solutions from solutions that will be performed by genetic operators
        index_aux = _get_int_random_number(&max_random);
        ind_1 = index_solutions_to_reproduce[index_aux];
        index_aux = _get_int_random_number(&max_random);
        ind_2 = index_solutions_to_reproduce[index_aux];        
        //Obtaining the proteins that were choosen from solutions
        prot_aux_1 = (protein_t*) solutions[ind_1].representation;
        prot_aux_2 = (protein_t*) solutions[ind_2].representation;        
        //Appling crossover operator between prot_aux_1 and prot_aux_2. Resulting pop_new[i]
        apply_crossover(&pop_new[i], prot_aux_1, prot_aux_2);        
        //Appling mutation operator in pop_new[i]        
        apply_mutation(&pop_new[i], in_para);
    }

}

int ea_mono(const input_parameters_t *in_para){

    char *prefix;
    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    
    protein_t *population_p; // main population
    protein_t *population_new; //population obtained by reprodution

    solution_t *solutions_p; // main solution
    int *index_solutions_to_reproduce; // stores the index of solutions that were selected to reproduce

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

    //Allocating PDB ATOMS
    population_p = allocateProtein(&in_para->size_population);
    population_new = allocateProtein(&in_para->size_population);
    
    //Loading initial population
    load_initial_population_file(population_p, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);    

    //Setting population_new
    copy_protein_population(population_new, population_p, &in_para->size_population);
    initialize_protein_population_atoms(population_new, &in_para->size_population);

    //Saving topology of population 
    prefix = Malloc(char,10);
    strcpy(prefix, "prot");
    save_topology_population(population_p, &in_para->size_population, 
    in_para->path_local_execute, prefix);    
    free(prefix);

/**************** STARTING Mono-Objetive Evolutionary Algorithm *************************/
    display_msg("Starting Mono-Objetive Evolutionary Algorithm \n");
    initialize_algorithm_execution(primary_sequence, in_para);
    init_gromacs_execution();
    index_solutions_to_reproduce = Malloc(int, in_para->size_population);
    //Setting solutions
    solutions_p = allocate_solution(&in_para->size_population, &in_para->number_fitness);        
    //Setting reference of proteins to solution 
    set_proteins2solutions(solutions_p, population_p, &in_para->size_population);
    //Computing objectives of solutions with GROMACS
    get_gromacs_objectives(solutions_p, in_para);
    for (int g = 1; g <= in_para->number_generation; g++){        
        //Seleting solutions (individuals) which will be reproduced. Now, all population will be reproduced
        select_individual_mono(index_solutions_to_reproduce, &in_para->size_population, 
            solutions_p, &in_para->size_population);
        //Reprodution of population - Building new proteins based on selected individuals
        reproduce_protein(population_new, solutions_p, index_solutions_to_reproduce, 
            in_para, &in_para->size_population);        
        //Updating the atoms of population_p from population_new
        copy_protein_population_atoms(population_p, population_new, &in_para->size_population);         
        /* Computing objectives of solutions with GROMACS since population_p was 
         * updated by new population (generated by reprodution) */
        get_gromacs_objectives(solutions_p, in_para);        
        //Saving information of generation
        update_execution_algorithms(solutions_p, &g);
    }
    // Sorting solutions_p
    // Saving sorted population
    finish_gromacs_execution();
    free(index_solutions_to_reproduce);
/**************** FINISHED Mono-objetive Evolutionary Algorithm *************************/     

    desallocate_solution(solutions_p, &in_para->size_population);
    desallocateProtein(population_p, &in_para->size_population);
    desallocate_primary_seq(primary_sequence);
    _finish_random_gsl();

	return 0;
}

