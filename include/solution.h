#ifndef OLD_SOLUTION_H
#define OLD_SOLUTION_H

#include "solution_types.h"

solution_t * allocate_solution(const int *nsol, const int *nbj);
void desallocate_solution(solution_t *sol, const int *nsol);
double get_specific_objective(const solution_t *sol, const int *obj);
double get_oposite_specific_objective(const solution_t *sol, const int *obj);

#endif