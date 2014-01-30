#include <stdlib.h>
#include <string.h>

#include "defines.h"
#include "enums.h"
#include "parameters_type.h"
#include "ea_mono.h"
#include "protein.h"
#include "topology.h"
#include "pdbio.h"
//#include "energy_gromacs.h"
#include "pdbatom.h"
#include "messages.h"
#include "ea_types.h"
#include "ea.h"
#include "string_owner.h"
#include "futil.h"
#include "random_number_gsl.h"

static void check_parameters(const input_parameters_t *in_para){
    /*Check the parameters to run EA mono*/
	if (in_para->number_individual_select_reproduce > in_para->size_population){
		fatal_error("The parameter NumberIndividualReproduce must be smaller than SizePopulation. Please check it!");
	}
}

int ea_mono(const input_parameters_t *in_para){
	check_parameters(in_para);
	//Popoulations of EA
	//protein **population_p, **population_new, **population_sel_repro;
    primary_seq_t *primary_sequence; // Primary Sequence of Protein

    fatal_error("Em topologia buscar os nomes dos atomos pelo seu nome e nao tipo o qual tenho de remover");

/*    
    display_msg("Build Topology \n");
	top_global = allocateTop_Global(&numatom, &nresiduos,&bond_angles,
			&num_protein_side_chains,&num_dihedral_angles_type);
    build_top_global(primary_sequence,top_global);
    save_topology(in_para->path_local_execute, in_para->top_file,top_global);

    display_msg("Allocating populations \n");
	population_p = allocatePopulation(in_para->size_population,nresiduos,
			in_para->number_fitness, top_global->numatom);
	population_new = allocatePopulation(in_para->size_population,nresiduos,
			in_para->number_fitness, top_global->numatom);
	population_sel_repro = allocatePopulation(in_para->number_individual_select_reproduce,
			nresiduos, in_para->number_fitness, top_global->numatom);

	display_msg("Reading population \n");
    readPopulationFile(population_p, &in_para->size_population,
    		in_para->path_local_execute,in_para->initial_pop_file_name,
    		primary_sequence,&nresiduos, top_global);

    display_msg("Starting Evolutionary Algorithm Mono-Objective \n");
    ea_initialize(primary_sequence, &nresiduos, top_global, in_para);
	initialize_fitness_in_file(in_para->path_local_execute,
			fitness_file_name);
    for (int g = 1; g <= in_para->number_generation;g++){
    	show_generation_msg(&g);
    	fitness_gromacs(population_p,in_para, top_global);
    	set_fitness_in_file(in_para->path_local_execute,fitness_file_name,
    			population_p,&in_para->size_population,&g);
    	select_individual_mono(population_sel_repro, population_p, &nresiduos,
    			in_para);
    	reproduce_population(population_new,population_sel_repro,in_para);
    	display_msg("Building random individuals \n");
    	build_random_individuals_mono(population_new,&in_para->size_population,
    			&in_para->number_individual_select_reproduce);
    	update_generation(population_p, &g);
    	/* population_p is updated until g < number_generation
    	 * because population_p is used to generate the final results and the
    	 * population_new will not be evaluated. So, We need to maintain the
    	 * population_p unchanged.
    	 */
/*    	 
    	if (g < in_para->number_generation){
    		copy_population_without_allocating(population_p,population_new,
    				&in_para->size_population);
    	}

    }
    display_msg("Generating Final Results \n");
    build_final_results(population_p,in_para, primary_sequence,
    		top_global);

    display_msg("Deallocating variables \n");
    ea_finished();
    deAllocateAmino(primary_sequence,nresiduos);
    desAllocateTop_Global(top_global);
    deAllocatePopulation(population_p, in_para->size_population);
    deAllocatePopulation(population_new, in_para->size_population);
    deAllocatePopulation(population_sel_repro, in_para->
    		number_individual_select_reproduce);
    finish_random();
*/    
	return 0;
}
