#include<stdlib.h>

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
