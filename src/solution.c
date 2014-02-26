#include <stdio.h>
#include <string.h>

#include "solution.h"
#include "defines.h"

/** Allocates an array of solution
* nsol is the number of solutions. Its number will be population size
* nbj is the number of objectives
*/
solution_t * allocate_solution(const int *nsol, const int *nbj){
	solution_t * solution_aux;	
	solution_aux = Malloc(solution_t,*nsol);
	//Each solution
	for (int s = 0; s < *nsol; s++){
		solution_aux[s].num_obj = *nbj;
		solution_aux[s].obj_values = Malloc(double, solution_aux[s].num_obj);
		solution_aux[s].representation = NULL;
	}
	return solution_aux;
}

void desallocate_solution(solution_t *sol, const int *nsol){
	for (int s = 0; s < *nsol; s++){
		free(sol[s].obj_values);
	}
	free(sol);
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