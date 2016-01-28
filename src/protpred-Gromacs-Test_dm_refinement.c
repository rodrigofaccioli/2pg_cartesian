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

int main(int argc, char *argv[])
{
	input_parameters_t *in_para;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);

	//MDP File
	char *mdp_test = Malloc(char, MAX_FILE_NAME);

	strcpy(mdp_test, "teste.mdp");
	

    primary_seq_t *primary_sequence = NULL; // Primary Sequence of Protein
    protein_t *population_p = NULL; //main population
    solution_t *solutions_p = NULL; //main solutions

    //Sample Test
    int one = 1;

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating population_p
    population_p = allocateProtein(&one);

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &in_para->size_population,
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);

    //Initialize Gromacs Execution
    init_gromacs_execution();

    //Test
    call_pdb2gmx_for_md(in_para->initial_pop_file_name, in_para->path_local_execute, in_para->path_gromacs_programs, in_para->force_field);
    call_mdrun2md(in_para->initial_pop_file_name, in_para->path_local_execute, in_para->path_gromacs_programs, mdp_test);
    call_trjconv2md(in_para->initial_pop_file_name, in_para->path_local_execute, in_para->path_gromacs_programs);

    //Finish Gromacs Execution
    finish_gromacs_execution();
    free(mdp_test);

    return 0;
}
