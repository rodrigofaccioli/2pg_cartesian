#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "random_number_gsl.h"
#include "defines.h"

long int _get_seed(){
	return get_seed();
}

int _get_int_random_number(__const int *max_number){
	if (max_number == NULL){
		int max_int = MAX_INT_NUMBER_RANDOM;
		return get_int_number(&max_int);
	}else{
		return get_int_number(max_number);
	}
}

double _get_double_random_number(){
	return get_double_number();
}

double _get_double_gama(const int *a, const int *b){
	return get_double_gama(a,b);
}

double _get_double_gauss(const int *g){
	return get_double_gauss(g);
}

void _finish_random_gsl(){
	finish_random();
}

/** Generates a psuedo-random float between 0.0 and 0.999 */
float _get_float(){ 
    return rand()/(float(RAND_MAX)+1); 
} 


/** Generates a psuedo-random float between 0.0 and max */
float _get_float_max(const float *max){
    return _get_float()* (*max); 
} 

/** Generates a psuedo-random float between min and max */
float _get_float_random_interval(const float *min, const float *max){
    if (*min > *max) 
        return _get_float()*(*min - *max) + *max;     
    else 
        return _get_float()*(*max - *min) + *min; 
}

