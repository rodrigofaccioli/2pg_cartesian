#ifndef OLD_EA_NSGA2_H
#define OLD_EA_NSGA2_H

#include "parameters_type.h"
#include "ea_nsga2_types.h"
#include "dominance_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
int ea_nsga2(const input_parameters_t *in_para);
ea_nsga2_t * allocate_nsga2(const input_parameters_t *in_para);
_2PG_CARTESIAN_EXPORT
ea_nsga2_t * allocate_nsga2_without_allocation_of_representation(const input_parameters_t *in_para);
_2PG_CARTESIAN_EXPORT
void desallocate_solution_nsga2(ea_nsga2_t *nsga2_sol, const int *size);
void set_nsga2_solution_in_solution(solution_t *solutions, 
    const ea_nsga2_t * nsga2_solutions, const int *size );
_2PG_CARTESIAN_EXPORT
void compute_fronts(ea_nsga2_t *nsga2_solutions, dominance_t * dominance,
        const int *size);    

#endif
