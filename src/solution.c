#include <stdio.h>
#include <string.h>

#include "solution.h"
#include "defines.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif


/** Initialize a Solution
* solution_aux represents the solution which will be initialized
* nbj represents the number of objectivies
*/
void initialize_solution(solution_t * solution_aux, const int *nbj){
		solution_aux->ID = -1;
		solution_aux->num_obj = *nbj;
		solution_aux->obj_values = Malloc(double, solution_aux->num_obj);
		solution_aux->representation = NULL;
}

/** Allocates an array of solution
* nsol is the number of solutions. Its number will be population size
* nbj is the number of objectives
*/
_2PG_CARTESIAN_EXPORT
solution_t * allocate_solution(const int *nsol, const int *nbj){
	solution_t * solution_aux;	
	solution_aux = Malloc(solution_t,*nsol);
	//Each solution
	for (int s = 0; s < *nsol; s++){
		initialize_solution(&solution_aux[s], nbj);
	}
	return solution_aux;
}

_2PG_CARTESIAN_EXPORT
void desallocate_solution(solution_t *sol, const int *nsol){
	for (int s = 0; s < *nsol; s++){
		free(sol[s].obj_values);
	}
	free(sol);
	sol = NULL;
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

/** Returns the displyed value of an specific objective */
double get_displayed_value_of_objective(const solution_t *sol, 
	const int *index, const int *obj, 
	const type_fitness_energies_t *fitness_energies){
	if ( (fitness_energies[*obj] == fit_hbond) ||
		(fitness_energies[*obj] == fit_hydrophilic)||
		(fitness_energies[*obj] == fit_hbond_main)	||
		(fitness_energies[*obj] == fit_stride_total) ||
		(fitness_energies[*obj] == fit_stride_helix)	||
		(fitness_energies[*obj] == fit_stride_beta) ) {
		return get_oposite_specific_objective(&sol[*index],obj);
	}else{
		return get_specific_objective(&sol[*index],obj);
	}
}

/** Returns index of solution based on value of objective
* sol contains all solutions
* size of sol
* obj_get_index index of objective that want to find out index of solution
* value is the value that want to find out index of solution
*/
int get_solution_index_by_objective_value(const solution_t *sol, const int *size,
            const int *obj_get_index, const double *value){
	for (int s = 0; s < *size; s++){
		if (sol[s].obj_values[*obj_get_index] == *value){
			return s;
		}
	}
	return -1;
}
