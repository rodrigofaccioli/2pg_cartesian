#include <stdlib.h>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "random_number_gsl.h"
#include "osutil.h"

extern unsigned long int gsl_rng_default_seed;

/*Must be started with NULL*/
static gsl_rng *r = NULL;
static unsigned long int s;
static const gsl_rng_type * T;

static void set_seed(){
	/* The seed value is based on process number */
	//gsl_rng_default_seed = get_pid();
	s = get_pid();//random();
	gsl_rng_env_setup();
	gsl_rng_set(r,s);
}

static void init_random(){
	/*Starts the random number generator.
	 * This module uses GNU Scientific Library (GSL)
	*/
	if (r == NULL){
		/*First execution of generator*/
		  gsl_rng_default_seed = get_pid()+random();
		  T = gsl_rng_default;
		  r = gsl_rng_alloc (T);
		//gsl_rng_env_setup();
		//r = gsl_rng_alloc(gsl_rng_ranlxs0); //gsl_rng_mt19937
	}
	//set_seed();
}

void finish_random()
{
       gsl_rng_free(r); 
}

long int get_seed(){
	/*Returns seed value*/
	return gsl_rng_default_seed;
}

int get_int_number(__const int *max){
	/*Returns a random integer number*/
	init_random();
	return gsl_rng_uniform_int(r,*max);
}

double get_double_number(){
	/*Returns a random double number*/
	init_random();
	return gsl_rng_uniform(r);
}


double get_double_gama(const int *a, const int *b){
	/*Returns a random double number based on Gama distribution*/
	init_random();
	return gsl_ran_gamma(r,*a,*b);
}

double get_double_gauss(const int *g){
	/*Returns a random double number based on Gauss distribution*/
	init_random();
	return gsl_ran_gaussian(r,*g);
}
