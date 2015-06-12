#include <stdlib.h>
#include <string.h>
#include <math.h>

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
#include "rotation.h"
#include "math_owner.h"
#include "vector_math.h"
#include "algorithms.h"

int main(int argc, char *argv[]){
	input_parameters_t *in_para;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);

    char *crossover_file_name = NULL;
    int ind, ind_ref_1, ind_ref_2;

    //Setting crossover file name 
    crossover_file_name = Malloc(char, MAX_FILE_NAME);
    strcpy(crossover_file_name, "pop_crossover.pdb");

    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    
    protein_t *population_p; // main population
    protein_t *new_population; // main population

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating PDB ATOMS of population_p
    population_p = allocateProtein(&in_para->size_population);
    new_population  = allocateProtein(&in_para->size_population);   

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);
    //Only to initialize new_population
    load_initial_population_file(new_population, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);
    initialize_protein_population_atoms(new_population, &in_para->size_population);

    for (ind = 0; ind < in_para->size_population; ind++){
        ind_ref_1 = 0;
        ind_ref_2 = 1;
        apply_crossover(&new_population[ind], &population_p[ind_ref_1], &population_p[ind_ref_2], in_para->crossovers);
    }

    save_population_file(new_population, in_para->path_local_execute, crossover_file_name, &in_para->size_population);

	deAllocateload_parameters(in_para);

	display_msg("Done 2PG Crossover \n");
	return 0;
}