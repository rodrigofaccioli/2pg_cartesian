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


int main(int argc, char *argv[]){
	input_parameters_t *in_para;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);
	
    primary_seq_t *primary_sequence = NULL; // Primary Sequence of Protein
    protein_t *population_p = NULL;

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating PDB ATOMS of population_p
    population_p = allocateProtein(&in_para->size_population);    

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);

    desallocateProtein(population_p, &in_para->size_population);
    deAllocateload_parameters(in_para);

	return 0;
}
