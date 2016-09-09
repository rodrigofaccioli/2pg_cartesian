#ifndef OLD_NSGA2_TYPE_H
#define OLD_NSGA2_TYPE_H

#include "dominancep_type.h"

typedef struct nsga2t{
	// Name of method's file or an ID that identifies a solution.
	char *file_name;
	// Front where the solution belongs to
	int front;
	// Crowding distance of the solution
	double crowding_distance;
	//Pointer to dominancep
	const dominancep_t *dominancep;
	// Ranking
	int ranking;
 }nsga2_t;


#endif
