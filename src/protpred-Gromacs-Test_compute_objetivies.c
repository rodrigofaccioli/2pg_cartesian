#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "load_parameters.h"
#include "messages.h"
#include "defines.h"
#include "solution.h"
#include "futil.h"
#include "string_owner.h"
#include "protein.h"
#include "topology.h"
#include "pdbio.h"
#include "pdbatom.h"
#include "aminoacids.h"
#include "aminoacids_io.h"
#include "populationio.h"
#include "topologyio.h"
#include "gromacs_objectives.h"
#include "objective.h"
#include "algorithms.h"

int main(int argc, char *argv[]){
	input_parameters_t *in_para;
    int ger = 0;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);

    int ind, obj;
    char *objective_name = NULL;
    objective_name = Malloc(char, MAX_RANDOM_STRING);

    primary_seq_t *primary_sequence = NULL; // Primary Sequence of Protein
    protein_t *population_p = NULL; //main population
    solution_t *solutions_p = NULL; //main solutions

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating population_p
    population_p = allocateProtein(&in_para->size_population);    

    //Allocating solutions
    solutions_p = allocate_solution(&in_para->size_population, &in_para->number_fitness);

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);

    //Setting reference of proteins to solution 
    set_proteins2solutions(solutions_p, population_p, &in_para->size_population);

    //Initialize algorithm execution
    initialize_algorithm_execution(primary_sequence, in_para);    

    //Initialize Gromacs Execution
    init_gromacs_execution();

    display_msg("Computing Objectives\n");
    get_gromacs_objectives(solutions_p, in_para);        

    //Finish Gromacs Execution
    finish_gromacs_execution();

    //Saving objetive files
    build_fitness_files(solutions_p, &ger, &in_para->size_population);

    free(objective_name);
    desallocate_solution(solutions_p, &in_para->size_population);
    desallocateProtein(population_p, &in_para->size_population);
    desallocate_primary_seq(primary_sequence);
    deAllocateload_parameters(in_para);

	return 0;
}
