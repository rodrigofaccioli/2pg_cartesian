#include <stdio.h>
#include <string.h>

#include "algorithms.h"
#include "algorithms_types.h"
#include "messages.h"
#include "defines.h"
#include "protein_type.h"
#include "populationio.h"
#include "solutionio.h"
#include "objective.h"


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

static void build_fitness_files(const solution_t *solutions, const int *generation,
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

