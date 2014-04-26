#ifndef OLD_DOMINANCE_TYPE_H
#define OLD_DOMINANCE_TYPE_H

#define MIN_DOMINATED 0
#define IS_DOMINANTED 0

#include "solution_types.h"

typedef struct sdominance{
	/* indicates the max index of set_dominance that are using. 
	 * It represents the number of solutions that are dominated 
	 * by solution (sol).
	*/
	int max_dominated;
	// Array of solutions that are dominated by solution (sol) 
	int *set_dominated;
	/* indicates how many solutions dominate solution (sol).
	 * When its value is zero means that sol is 
	 * not dominated. Therefore, it is a good solution. 
	*/
	int how_many_solutions_dominate_it;
	//Pointer to solution
	const solution_t *sol;
 }dominance_t;


#endif