#ifndef OLD_DOMINANCEP_TYPE_H
#define OLD_DOMINANCEP_TYPE_H

#include "solution_types.h"

typedef struct sdominancep{
	// Amount of solutions dominated by solution (sol)
	int amount_dominated;
	// Array of solutions that are dominated by solution (sol) 
	int *set_dominated;
	// Amount of solutions that dominate solution (sol)
	int amount_dominates_me;
	// Array of solutions that dominate solution (sol) 
	int *set_dominates_me;
	//Pointer to solution (sol)
	const solution_t *sol;
 }dominancep_t;


#endif
