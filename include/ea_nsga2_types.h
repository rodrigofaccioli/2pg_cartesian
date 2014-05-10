#ifndef OLD_EA_NSGA2_TYPES_H
#define OLD_EA_NSGA2_TYPES_H

#include "solution_types.h"

typedef struct sea_nsga2{	
	int front; //Represents front of solution (sol)
	float crowding_distance; //Represents crowding distance of solution (sol)
	solution_t *sol; //Represents a solution.
}ea_nsga2_t;


#endif
