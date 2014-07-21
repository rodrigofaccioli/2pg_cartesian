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
#include "objective.h"
#include "solutionio.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

static int compare_solution(const void *x, const void *y){
    double fx, fy;
    fx = ((solution_t *)x)->obj_values[0];
    fy = ((solution_t *)y)->obj_values[0];
    if (fx > fy){
        return 1;
    }else if ( (fx -fy) == 0.0001){
        return 0;
    }else{
        return -1;
    }
}

static void build_final_results(solution_t *pop,const input_parameters_t *in_para ){

    char *sorted_pop_file, *sorted_fit_file;
    int obj;
    protein_t *population_aux;
    const protein_t *prot;


    sorted_pop_file = Malloc(char, MAX_FILE_NAME);
    sorted_fit_file = Malloc(char, MAX_FILE_NAME);
    population_aux = allocateProtein(&in_para->size_population);

    obj = 0;
    strcpy(sorted_pop_file, "pop_sorted_file.pdb");
    type_fitness_energies2str(sorted_fit_file, &in_para->fitness_energies[obj]);
    strcat(sorted_fit_file,"_sort.fit");

    //sorting population
    qsort (pop, in_para->size_population, sizeof (solution_t), compare_solution);
    //Saving population
    for (int p = 0; p < in_para->size_population; p++){
        prot = (protein_t*) pop[p].representation;
        population_aux[p].p_atoms = allocate_pdbatom(&prot->p_topol->numatom);
        copy_protein(&population_aux[p], prot);
    }
    save_population_file(population_aux, in_para->path_local_execute, sorted_pop_file, 
        &in_para->size_population);
    //Saving objective values
    save_solution_file(in_para->path_local_execute, sorted_fit_file, &obj, 
            pop, &in_para->size_population, &in_para->number_generation, in_para);

    for (int i = 0; i < in_para->size_population; i++){
        desAllocate_pdbatom(population_aux[i].p_atoms);
    }   
    desallocateProtein(population_aux, &in_para->size_population);
    free(sorted_pop_file); 
    free(sorted_fit_file);
}


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
        apply_crossover(&pop_new[i], prot_aux_1, prot_aux_2, in_para->crossovers);
        //Appling mutation operator in pop_new[i]        
        apply_mutation(&pop_new[i], in_para);
    }

}

_2PG_CARTESIAN_EXPORT
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
    for (int i = 0;  i < in_para->size_population; i++){
        population_new[i].p_atoms = allocate_pdbatom(&population_p[0].p_topol->numatom);            
    }
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
    build_final_results(solutions_p, in_para);    
    finish_gromacs_execution();
    free(index_solutions_to_reproduce);
/**************** FINISHED Mono-objetive Evolutionary Algorithm *************************/     
    
    desallocate_solution(solutions_p, &in_para->size_population);
    desallocateProtein(population_p, &in_para->size_population);
    desallocate_primary_seq(primary_sequence);

    return 0;    
}

