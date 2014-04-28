#include <stdio.h>
#include <string.h>

#include "solution.h"
#include "defines.h"


/** Initialize a Solution
* solution_aux represents the solution which will be initialized
* nbj represents the number of objectivies
*/
void initialize_solution(solution_t * solution_aux, const int *nbj){
		solution_aux->num_obj = *nbj;
		solution_aux->obj_values = Malloc(double, solution_aux->num_obj);
		solution_aux->representation = NULL;
}

/** Allocates an array of solution
* nsol is the number of solutions. Its number will be population size
* nbj is the number of objectives
*/
solution_t * allocate_solution(const int *nsol, const int *nbj){
	solution_t * solution_aux;	
	solution_aux = Malloc(solution_t,*nsol);
	//Each solution
	for (int s = 0; s < *nsol; s++){
		initialize_solution(&solution_aux[s], nbj);
	}
	return solution_aux;
}

void desallocate_solution(solution_t *sol, const int *nsol){
	for (int s = 0; s < *nsol; s++){
		free(sol[s].obj_values);
	}
	free(sol);
}

/** Copy all objectivies between two solutions
*/ 
void copy_solution_objectivies(solution_t *dest, const solution_t *source){
    dest->num_obj = source->num_obj;
    for (int obj = 0; obj < dest->num_obj; obj++){
    	dest->obj_values[obj] = source->obj_values[obj];
    }
}

/** Returns the value of an specific objective */
double get_specific_objective(const solution_t *sol, const int *obj){
	return sol->obj_values[*obj];
}


/** Returns the oposite value of an specific objective */
double get_oposite_specific_objective(const solution_t *sol, const int *obj){
	/* objectives that must be maximized when they are obtained these values
	 * are multiplied by -1 because 2PG works considerating minimization.
	 * However, when they will be stored, they must be written in original
	 * value 
	 */	
	return sol->obj_values[*obj] * (-1);
}