#ifndef OLD_SPEA2_TYPE_H
#define OLD_SPEA2_TYPE_H

#include "dominancep_type.h"

typedef struct spea2t{
	// Name of method's file or an ID that identifies a solution.
	char *file_name;
	// Raw fitness value as in SPEA2 
	double raw_fitness;
	// Density value as in SPEA2
	double density;
	// Fitness value as the sum of raw fitness and density just like SPEA2
	double fitness;
	//Pointer to dominancep
	const dominancep_t *dominancep;
	// Ranking
	int ranking;
 }spea2_t;


#endif
