#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "defines.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
int _get_int_random_number(__const int *max_number){
	int max_int;
	srand(time(NULL));
	if (max_number == NULL){
		max_int = MAX_INT_NUMBER_RANDOM;		
	}else{
		max_int = *max_number;		
	}
	return rand() % max_int;
}

/** Generates a psuedo-random float between 0.0 and 0.999 */
_2PG_CARTESIAN_EXPORT
float _get_float(){ 
    return rand()/(float(RAND_MAX)+1); 
} 

/** Generates a psuedo-random float between 0.0 and max */
_2PG_CARTESIAN_EXPORT
float _get_float_max(const float *max){
    return _get_float()* (*max); 
} 

/** Generates a psuedo-random float between min and max */
_2PG_CARTESIAN_EXPORT
float _get_float_random_interval(const float *min, const float *max){
    if (*min > *max) 
        return _get_float()*(*min - *max) + *max;     
    else 
        return _get_float()*(*max - *min) + *min; 
}
