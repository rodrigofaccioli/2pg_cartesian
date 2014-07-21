#ifndef OLD_ALGORITHMS_H
#define OLD_ALGORITHMS_H

#include "parameters_type.h"
#include "aminoacids_types.h"
#include "solution_types.h"
#include "protein_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
void initialize_algorithm_execution(const primary_seq_t *primary_sequence_aux,
		const input_parameters_t *in_para_aux);
_2PG_CARTESIAN_EXPORT
void update_execution_algorithms(const solution_t *solutions, const int *tag);
_2PG_CARTESIAN_EXPORT
int get_choose_residue(const int *num_res_prot);
void apply_crossover(protein_t *ind_new, const protein_t *prot_1, 
    const protein_t * prot_2, type_crossoers_t *crossovers);
void apply_mutation(protein_t *ind_new, const input_parameters_t *in_para);
void build_fitness_files(const solution_t *solutions, const int *generation,
        const int *pop_size);
int get_started_generation(const int *start);
#endif
