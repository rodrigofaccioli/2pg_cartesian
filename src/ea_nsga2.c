#include <stdlib.h>
#include <string.h>


#include "ea_nsga2.h"
#include "messages.h"
#include "defines.h"
#include "enums.h"
#include "protein.h"
#include "topology.h"
#include "pdbio.h"
#include "pdbatom.h"
#include "messages.h"
#include "algorithms.h"
#include "string_owner.h"
#include "futil.h"
#include "aminoacids.h"
#include "aminoacids_io.h"
#include "populationio.h"
#include "topology.h"
#include "topologyio.h"
#include "rotation.h"
#include "math_owner.h"
#include "solution.h"
#include "gromacs_objectives.h"
#include "algorithms.h"
#include "randomlib.h"
#include "objective.h"
#include "solutionio.h"
#include "dominance.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

#define INIT_FRONT -1

/* It is used in compare_objective function which is used in qsort
* It is set in set_crowding_distance function.
*/
static int index_obj_compare_objective = -1;

/* It is used in compare_objective_solution function which is used in qsort
* It is used for sorting solution by first objective
*/
static int index_obj_compare_objective_solution = -1;

/* It represents global topology that is used in 
* nsga2 solutions
*/
static top_global_t *nsga2_global_topol = NULL;
static protein_t *protein_ref_global = NULL;

ea_nsga2_t * allocate_nsga2(const input_parameters_t *in_para){
	ea_nsga2_t * aux = NULL;
    protein_t *prot_aux = NULL;
    int prot = 1;
	aux = Malloc(ea_nsga2_t,in_para->size_population);
	for (int ind = 0; ind < in_para->size_population; ind++){
		aux[ind].front = INIT_FRONT;
		aux[ind].crowding_distance = MIN_DIST;		
		aux[ind].sol = Malloc(solution_t,1);
		initialize_solution(aux[ind].sol, &in_para->number_fitness);
        prot_aux = allocateProtein(&prot);
        prot_aux->p_atoms = allocate_pdbatom(&nsga2_global_topol->numatom);
        prot_aux->p_topol = nsga2_global_topol;
        aux[ind].sol->representation = prot_aux;

	}
	return aux;
}

_2PG_CARTESIAN_EXPORT
ea_nsga2_t * allocate_nsga2_without_allocation_of_representation(const input_parameters_t *in_para){
    ea_nsga2_t * aux = NULL;
    protein_t *prot_aux = NULL;
    int prot = 1;
    aux = Malloc(ea_nsga2_t,in_para->size_population);
    for (int ind = 0; ind < in_para->size_population; ind++){
        aux[ind].front = INIT_FRONT;
        aux[ind].crowding_distance = MIN_DIST;      
        aux[ind].sol = Malloc(solution_t,1);
        initialize_solution(aux[ind].sol, &in_para->number_fitness);        
    }
    return aux;
}


ea_nsga2_t * allocate_nsga2_RT(const int *size, const input_parameters_t *in_para){    
    ea_nsga2_t * aux;
    protein_t *prot_aux = NULL;
    int prot = 1;    
    aux = Malloc(ea_nsga2_t, *size);
    for (int ind = 0; ind < *size; ind++){
        aux[ind].front = INIT_FRONT;
        aux[ind].crowding_distance = MIN_DIST;      
        aux[ind].sol = Malloc(solution_t,1);
        initialize_solution(aux[ind].sol, &in_para->number_fitness);
        prot_aux = allocateProtein(&prot);
        prot_aux->p_atoms = allocate_pdbatom(&nsga2_global_topol->numatom);
        prot_aux->p_topol = nsga2_global_topol;
        aux[ind].sol->representation = prot_aux;

    }
    return aux;
}

_2PG_CARTESIAN_EXPORT
void desallocate_solution_nsga2(ea_nsga2_t *nsga2_sol, const int *size){
    int size_aux = 1;    
    protein_t *prot_aux;
    for (int i = 0; i < *size; i++){
        if (nsga2_sol[i].sol->representation != NULL){
            prot_aux = (protein_t*) nsga2_sol[i].sol->representation;
            desAllocate_pdbatom(prot_aux->p_atoms);
            prot_aux->p_topol = NULL;  
            desallocateProtein(prot_aux, &size_aux);
            nsga2_sol[i].sol->representation = NULL;
        }        
        desallocate_solution(nsga2_sol[i].sol, &size_aux);
    }
    free(nsga2_sol);
    nsga2_sol = NULL;
}


void desallocate_solution_nsga2_RT(ea_nsga2_t *nsga2_sol, const int *size){
    int size_aux = 1;    
    protein_t *prot_aux;
    for (int i = 0; i < *size; i++){
        if (nsga2_sol[i].sol->representation != NULL){
            prot_aux = (protein_t*) nsga2_sol[i].sol->representation;   
            desAllocate_pdbatom(prot_aux->p_atoms);
            prot_aux->p_topol = NULL;  
            desallocateProtein(prot_aux, &size_aux);
            nsga2_sol[i].sol->representation = NULL;
        }
        desallocate_solution(nsga2_sol[i].sol, &size_aux);
    }
    free(nsga2_sol);
    nsga2_sol = NULL;    
}

void show_nsga2_solutions(const ea_nsga2_t *nsga2_sol, const int *size){
    for (int i = 0; i < *size; i++){
        printf("Index %i \n", i);
        printf("front %i \n", nsga2_sol[i].front); 
    }
}

/** Save PDB files from nsga2 solutions
* solutions represents solutions that want to save
* path where of will save
* file_name name of file
* size number of solutions
*/
void save_nsga2_population_file(const ea_nsga2_t *solutions, const char *path, 
    const char *file_name, const int *size){

    FILE * pdbfile; 
    const protein_t *p_curr;
    char *fname = path_join_file(path, file_name);    
    for (int s = 1; s <= *size; s++){
        p_curr = (protein_t*) solutions[s-1].sol->representation;
        if (s == 1){
            pdbfile = open_file(fname, fWRITE);
            writeHeader(pdbfile, 0.00, &p_curr->p_topol->numatom);
        }else{
            pdbfile = open_file(fname, fAPPEND);
        }
        writeModel(pdbfile, &s);
        writeATOM(pdbfile, p_curr->p_atoms, &p_curr->p_topol->numatom);
        writeEndModel(pdbfile);
        fclose(pdbfile);
    }
    free(fname);
}

/** Save PDB files from solution_t solutions
* solutions represents solutions that want to save
* path where of will save
* file_name name of file
* size number of solutions
*/
void save_solution_population_file(solution_t *solutions, const char *path, 
    const char *file_name, const int *size){

    FILE * pdbfile; 
    const protein_t *p_curr;
    char *fname = path_join_file(path, file_name);    
    for (int s = 1; s <= *size; s++){
        p_curr = (protein_t*) solutions[s-1].representation;
        if (s == 1){
            pdbfile = open_file(fname, fWRITE);
            writeHeader(pdbfile, 0.00, &p_curr->p_topol->numatom);
        }else{
            pdbfile = open_file(fname, fAPPEND);
        }
        writeModel(pdbfile, &s);
        writeATOM(pdbfile, p_curr->p_atoms, &p_curr->p_topol->numatom);
        writeEndModel(pdbfile);
        fclose(pdbfile);
    }
    free(fname);
}


/** Setting the protein representation in ea_nsga2_t struct
* nsga2_sol represents the ea_nsga2_t 
* pop population of protein
* pop_size size of population
*/
void set_proteins2nsga2_solutions(ea_nsga2_t *nsga2_sol, protein_t *pop, 
	const int *pop_size){
    protein_t *prot_aux_dest;
	for (int i = 0; i < *pop_size; i++){
        prot_aux_dest = (protein_t*) nsga2_sol[i].sol->representation;
		copy_protein(prot_aux_dest, &pop[i]);
	}
}

/** Copy the protein representation from ea_nsga2_t to solution
* p_solp is a solution population 
* nsga2_sol represents the ea_nsga2_t 
* pop_size size of population
*/
void copy_nsga2_solutions2solution(solution_t *p_sol,
    const ea_nsga2_t *nsga2_sol,  const int *size){

    const protein_t *prot_aux_source;
    protein_t *prot_aux_dest;
    int one = 1;    

    for (int i = 0; i < *size; i++){
        if  (p_sol[i].representation == NULL){
            prot_aux_dest = allocateProtein(&one);
            prot_aux_dest->p_atoms = allocate_pdbatom(&nsga2_global_topol->numatom);
            prot_aux_dest->p_topol = nsga2_global_topol;
            p_sol[i].representation = prot_aux_dest;            
        }
        //Copy protein
        prot_aux_source = (protein_t*) nsga2_sol[i].sol->representation;
        prot_aux_dest = (protein_t*) p_sol[i].representation;        
        copy_protein(prot_aux_dest, prot_aux_source);
    }
}


/** Computes each objective by GROMACS
* nsga2_solutions represents NSGA2 solutions
* in_para input parameters
*/
void nsga2_compute_objectives_with_gromacs(ea_nsga2_t *nsga2_solutions,
	const int *size, const input_parameters_t *in_para){

	for (int ind = 0; ind < *size; ind++){        
		get_gromacs_objectives_of_solution(nsga2_solutions[ind].sol, 
			in_para, &ind);
	}
}

/** Setting  nsga2_solutions in solution. 
* solutions is the solutions that will receive the objectivies
* nsga2_solutions contains the solutions that have values of objectivies
* size number of solutions
*/
void set_nsga2_solution_in_solution(solution_t *solutions, 
    const ea_nsga2_t * nsga2_solutions, const int *size ){

    for (int s = 0; s < *size; s++){        
        //Copy objectivies from nsga2_solutions to solutions
        for (int obj =0; obj < solutions[s].num_obj; obj++){
            solutions[s].obj_values[obj] = nsga2_solutions[s].sol->obj_values[obj];
        }
    }

}


/** Looking for front equal -1. It represents which is necessary to set front
* nsga2_solutions solutions that will set front
* size number of solutions 
*/
static boolean_t continue_front_computation(const ea_nsga2_t *nsga2_solutions,
        const int *size){   
    boolean_t ret = bfalse;
    for (int i = 0; i < *size; i++){
        if (nsga2_solutions[i].front == INIT_FRONT){
            ret = btrue;
            break;
        }
    }
    return ret;
}


/**
*/
static int get_how_many_zero_front_init(const ea_nsga2_t *nsga2_solutions, const dominance_t * dominance, const int *size){
    int ret;
    for (int i = 0; i < *size; i++){
        if ((dominance[i].how_many_solutions_dominate_it == 0) &&
            (nsga2_solutions[i].front == INIT_FRONT) ){
            ret = ret + 1;
        }        
    }
    return ret;
}

/** Computes fronts to nsga2 solutions based on dominance criterion
* nsga2_solutions who will receive computed front
* dominance who has computed the dominance criterion 
* size is number of solutions
*/
_2PG_CARTESIAN_EXPORT
void compute_fronts(ea_nsga2_t *nsga2_solutions, dominance_t * dominance,
        const int *size){    

    int front, min_dominates;
    int index_dominated_solutions;
    int index;
    int how_many_is_zero;
    int *temp_aux=NULL;

    //Initialize variables
    front = 0;
    min_dominates = 0;
    index_dominated_solutions = -1;    
    how_many_is_zero = -1;

    //Initialize fronts
    for (int i = 0; i < *size; i++){
        nsga2_solutions[i].front = INIT_FRONT;
    }

    //Check computation of fronts. While front equal -1, computes fronts
    while (continue_front_computation(nsga2_solutions, size) == btrue){
        how_many_is_zero = get_how_many_zero_front_init(nsga2_solutions, dominance, size) - 1;
        temp_aux = Malloc(int, how_many_is_zero);
        //Setting index for temp_aux based on nsga2_solutions
        int j = 0;
        for (int i = 0; i < *size; i++){ //Loop all popupation looking for index
            if ((dominance[i].how_many_solutions_dominate_it == 0) &&
                (nsga2_solutions[i].front == INIT_FRONT) ){
                temp_aux[j] = i;
                j = j + 1;
            }
        }
        //Setting front based on temp_aux
        for (int z = 0; z < how_many_is_zero; z++){
            index = temp_aux[z]; //Assign index obtained in before loop
            nsga2_solutions[index].front = front;
            //Update how_many_solutions_dominate_it of solutions that are dominated by i
            for (int d = 0; d < dominance[index].max_dominated; d++){
                index_dominated_solutions = dominance[index].set_dominated[d]; //Set of Solutions that is dominated by i
                if  ( (index_dominated_solutions >= 0) &&
                        (index_dominated_solutions < *size) ){
                    dominance[index_dominated_solutions].how_many_solutions_dominate_it = dominance[index_dominated_solutions].how_many_solutions_dominate_it -1;
                }else{
                    fatal_error("Index index_dominated_solutions is wrong \n");
                }
            }            
        }
        free(temp_aux);
        front = front +1;
    }
}

/** Setting front to nsga2_solutions_p
* nsga2_solutions_p nsga_2 solution_p
* size number of solutions
*/
void set_front_solution_p(ea_nsga2_t *nsga2_solutions_p, const int *size){
	dominance_t * dominance;
    solution_t *solutions;
	int num_obj;

    /*Based on first solution is got the number of objective since 
    * it is the same for all solutions
    */
    num_obj = nsga2_solutions_p[0].sol->num_obj;    

	dominance = allocate_dominance(size);
    solutions = allocate_solution(size, &num_obj);

    //Coping values of objective    
    set_nsga2_solution_in_solution(solutions, nsga2_solutions_p, size );
    //Setting dominance
	set_dominance(dominance, solutions, size);    
    //Setting front based on dominance concept
    compute_fronts(nsga2_solutions_p, dominance, size);

	desallocate_dominance(dominance, size);
    desallocate_solution(solutions, size);

}

const ea_nsga2_t * tournament_front(const ea_nsga2_t *s1, const ea_nsga2_t *s2){
    if (s1->front < s2->front){
        return s1;
    }else{
        return s2;
    }
}

/** Applies the reprodution in proteins that are into solutions
* pop_new is the new population. It will be generated by genetic operators
* solutions represent the current solutions. In these solutions is applied genetic operators 
* in_para input parameters
* Important: When crossover is set to none, new individual is obtained by coping one
* individual. Mutation operator is always employed.
*/
void reproduce_protein(protein_t *pop_new, const ea_nsga2_t *solutions,  
    const input_parameters_t *in_para){

    const ea_nsga2_t *nsga2_sol_1, *nsga2_sol_2; 
    protein_t *prot_aux_1, *prot_aux_2;
    int max_random, index_ind_1, index_ind_2;

    max_random = in_para->size_population -1;
    for (int i = 0; i < in_para->size_population; i++){
        //Getting first individual
        index_ind_1 = _get_int_random_number(&max_random);
        index_ind_2 = _get_int_random_number(&max_random);
        nsga2_sol_1 = tournament_front(&solutions[index_ind_1], 
            &solutions[index_ind_2]);
        prot_aux_1 = (protein_t*) nsga2_sol_1->sol->representation;
        //Checking if apply crossover
        if (in_para->crossovers[0] != crossoer_none){
            //Getting second individual
            index_ind_1 = _get_int_random_number(&max_random);
            index_ind_2 = _get_int_random_number(&max_random);
            nsga2_sol_2 = tournament_front(&solutions[index_ind_1], 
                &solutions[index_ind_2]);
            prot_aux_2 = (protein_t*) nsga2_sol_2->sol->representation;        
            // Appling crossover operator between prot_aux_1 and prot_aux_2.             
            apply_crossover(&pop_new[i], prot_aux_1, prot_aux_2, in_para->crossovers);
        }else{
            //Coping prot_aux_1 to pop_new[i] when is not used crossover. 
            copy_protein(&pop_new[i], prot_aux_1);
        }        
        //Appling mutation operator in pop_new[i]
        apply_mutation(&pop_new[i], in_para);        
    }
}

/** Copy all information between two nsga2 solutions
*/
void copy_nsga2_solution(ea_nsga2_t * dest, const ea_nsga2_t * source){
    const protein_t *prot_aux_source;
    protein_t *prot_aux_dest;

    dest->front             = source->front;
    dest->crowding_distance = source->crowding_distance;
    copy_solution_objectivies(dest->sol, source->sol);
    //Copy protein
    prot_aux_source = (protein_t*) source->sol->representation;
    prot_aux_dest   = (protein_t*) dest->sol->representation;
    copy_protein(prot_aux_dest, prot_aux_source);
}

void set_soluton_rt(ea_nsga2_t *solutions_rt, const ea_nsga2_t *solutions_p, 
    const ea_nsga2_t *solutions_q, const int *size){
    int index_rt = -1;
    //Copy solutions_p to solutions_rt
    for (int i = 0; i < *size; i++){
        index_rt++;
        copy_nsga2_solution(&solutions_rt[index_rt], &solutions_p[i]);
    }

    //Copy solutions_q to solutions_rt
    for (int i = 0; i < *size; i++){
        index_rt++;
        copy_nsga2_solution(&solutions_rt[index_rt], &solutions_q[i]);
    }
}

/** Compare two nsga2 solutions. It is used in
* set_dominance_and_crowding_distance_in_soluton_rt function
* to sort the solutions by front. It is q_sort command. 
*/
static int compare_front(const void *x, const void *y){
    int fx, fy;
    fx = ((ea_nsga2_t *)x)->front;
    fy = ((ea_nsga2_t *)y)->front;
    if (fx > fy){
        return 1;
    }else {
        return 0;
    }
}

/** Computes how many individuals there is in front 
*/
int compute_how_many_front(const ea_nsga2_t * solutions, 
            const int *size, const int *front_ref){
    int num = 0;
    for (int i = 0; i < *size; i++){
        if (solutions[i].front == *front_ref){
            num = num + 1;
        }
    }
    return num;
}


static int compare_objective_solution( const void *x, const void *y){
    solution_t * fx = (solution_t*) x;
    solution_t * fy = (solution_t*) y;

    if (fx->obj_values[index_obj_compare_objective_solution] < fy->obj_values[index_obj_compare_objective_solution]){
        return -1;
    }else if (fx->obj_values[index_obj_compare_objective_solution] == fy->obj_values[index_obj_compare_objective_solution]){
        return 0;
    }else if (fx->obj_values[index_obj_compare_objective_solution] > fy->obj_values[index_obj_compare_objective_solution]){
        return 1;    
    }
}

static int compare_objective( const void *x, const void *y){
    ea_nsga2_t * fx = (ea_nsga2_t*) x;
    ea_nsga2_t * fy = (ea_nsga2_t*) y;

    if (fx->sol->obj_values[index_obj_compare_objective] > fy->sol->obj_values[index_obj_compare_objective]){
        return 1;
    }else {
        return 0;
    }
}

void set_crowding_distance(ea_nsga2_t *solution, const int *size,
        const int *num_obj){

    //Initialize crowding distance
    for (int i = 0; i < *size; i++){
        solution[i].crowding_distance = 0.0;
    }

    for (int i_obj = 0; i_obj < *num_obj; i_obj++){
        index_obj_compare_objective = i_obj;
        //Sort population based on objective
        qsort (solution, *size , sizeof (ea_nsga2_t),
                      compare_objective);
        //Set the limits. It represents infinity
        solution[0].crowding_distance = MAX_DIST;
        solution[*size -1].crowding_distance = MAX_DIST;
        for (int i = 1; i < *size-1; i++){
            solution[i].crowding_distance += solution[i+1].sol->obj_values[index_obj_compare_objective] - solution[i-1].sol->obj_values[index_obj_compare_objective];
        }
    }
}


void compute_crowding_distance(ea_nsga2_t * solutions_rt, 
            const int *size, const input_parameters_t *in_para){

    int front, j;
    int how_many_ind_front;
    int index_p,index_rt;
    int index_begin_back_to_rt;
    int total_ind;

    ea_nsga2_t *temp_pop;

    front = 0;
    index_begin_back_to_rt = 0;
    total_ind = 0;

    how_many_ind_front = compute_how_many_front(solutions_rt, size, &front);
    while (how_many_ind_front > 0){
        //Compute the number of individuals that will use
        total_ind = index_begin_back_to_rt + how_many_ind_front;
        //Copy Individuals of front to temp_pop
        temp_pop = allocate_nsga2_RT(&how_many_ind_front, in_para);
        j = 0;
        //Obtain the individuals by front
        for (int index_rt = index_begin_back_to_rt; index_rt <total_ind; index_rt++){
            if (solutions_rt[index_rt].front == front){
                copy_nsga2_solution(&temp_pop[j], &solutions_rt[index_rt]);
                j++;
            }else{
                fatal_error("In compute_crowding_distance function the population must be sorted by front \n");
            }
        }

        //Apply to Crowding distance in temp_pop
        set_crowding_distance(temp_pop, &how_many_ind_front, &in_para->number_fitness);

        //Copy to solutions_rt from temp_pop
        j = 0;
        for (int index_rt = index_begin_back_to_rt; index_rt <total_ind; index_rt++){
            copy_nsga2_solution(&solutions_rt[index_rt], &temp_pop[j]);
            j++;
        }

        //Free temp_pop
        desallocate_solution_nsga2_RT(temp_pop,&how_many_ind_front);

        //Updates index_begin_back_to_rt
        index_begin_back_to_rt = index_begin_back_to_rt + how_many_ind_front;

        //Updates front
        front = front +1;

        // Computes how many individuals has for next front
        how_many_ind_front = compute_how_many_front(solutions_rt, size, &front);
    }
}

void set_dominance_and_crowding_distance_in_soluton_rt(ea_nsga2_t * solutions_rt, 
            const int *size_RT, const input_parameters_t *in_para){
    dominance_t * dominance;
    solution_t *solutions;
    int num_obj;

    /*Based on first solution is got the number of objective since 
    * it is the same for all solutions
    */
    num_obj = solutions_rt[0].sol->num_obj;    

    dominance = allocate_dominance(size_RT);
    solutions = allocate_solution(size_RT, &num_obj);

    //Coping values of objective    
    set_nsga2_solution_in_solution(solutions, solutions_rt, size_RT);
    //Setting dominance
    set_dominance(dominance, solutions, size_RT);    
    //Setting front based on dominance concept
    compute_fronts(solutions_rt, dominance, size_RT);
    //Sorting by front
    qsort(solutions_rt, *size_RT,  sizeof (ea_nsga2_t), compare_front);
    //Computes Crowding Distance
    compute_crowding_distance(solutions_rt, size_RT, in_para);

    desallocate_dominance(dominance, size_RT);
    desallocate_solution(solutions, size_RT);

}

static void set_objective_file_name_non_dominated(char *objective_file_name_non_dominated,
    const int *fit, const int *generation, const input_parameters_t *in_para){
    char fitness_name[MAX_RANDOM_STRING];
    char sger[MAX_RANDOM_STRING];
    type_fitness_energies2str(fitness_name, &in_para->fitness_energies[*fit]);
    sprintf(sger, "%d", *generation);
    strcat(fitness_name,"_");
    strcat(fitness_name,"NON_DOMINATED_");
    strcat(fitness_name,sger);
    sprintf(objective_file_name_non_dominated,"%s.fit",fitness_name );
}

static void build_objective_files_non_dominated(const solution_t *solutions, 
    const int *generation, const int *pop_size, const int *num_obj, 
    const input_parameters_t *in_para){
    int f;
    char *objective_file_name_non_dominated=NULL;
    for (f =0; f < *num_obj; f++){
        objective_file_name_non_dominated = Malloc(char, MAX_RANDOM_STRING);
        set_objective_file_name_non_dominated(objective_file_name_non_dominated, 
            &f,generation, in_para);
        save_solution_file(in_para->path_local_execute, 
            objective_file_name_non_dominated, &f, 
            solutions, pop_size, generation, in_para);
        free(objective_file_name_non_dominated);
    }
}

/** Returns index of nsga2 solution based on value of objective
* sol contains all solutions
* size of sol
* obj_get_index index of objective that want to find out index of solution
* value is the value that want to find out index of solution
*/
int get_nsga2_solution_index_by_objective_value(const ea_nsga2_t *sol, const int *size,
            const int *obj_get_index, const double *value){
    for (int s = 0; s < *size; s++){
        if (sol[s].sol->obj_values[*obj_get_index] == *value){
            return s;
        }
    }
    return -1;
}

void build_plot_xvg_file_initial_population(const solution_t *solutions_sorted, 
    const int *size, const int *num_obj, const char *local_execute, 
    const int *front, const int *ger, 
    const type_fitness_energies_t *fitness_energies, 
    const int *obj_get_index){

    char *file_name;
    char *xvg_line, *aux_str;
    double v_aux;
    int ind_id_before_sorted;

    file_name = Malloc(char, MAX_FILE_NAME);
    xvg_line = Malloc(char, MAX_LINE_FILE);
    aux_str = Malloc(char, 60);
    sprintf(file_name,"plot_nfront_%d_%d.xvg",*front, *ger);
    char *fname = path_join_file(local_execute,file_name);
    FILE *xvg_file = open_file(fname,fWRITE);
    fprintf (xvg_file,"#Obj1\tObj2\t\n");
    for (int i=0; i < *size; i++){
        for (int obj = 0; obj < *num_obj; obj++){
            v_aux = get_displayed_value_of_objective(solutions_sorted, &i, &obj, 
                fitness_energies);            
            sprintf(aux_str,"%f",v_aux);
            if (obj == 0){
                strcpy(xvg_line, aux_str);
            }else{
                strcat(xvg_line, aux_str);
            }
            strcat(xvg_line, "\t");
        }
        fprintf (xvg_file,"%s\n", xvg_line);
    }
    fclose(xvg_file);

    free(aux_str);
    free(xvg_line);
    free(file_name);
    free(fname);
}


void build_plot_xvg_file(const solution_t *solutions_sorted, 
    const int *size, const int *num_obj, const char *local_execute, 
    const ea_nsga2_t *solutions_full, const int *size_full, 
    const int *front, const int *ger, 
    const type_fitness_energies_t *fitness_energies, 
    const int *obj_get_index){

    char *file_name;
    char *xvg_line, *aux_str;
    double v_aux;
    int ind_id_before_sorted;

    file_name = Malloc(char, MAX_FILE_NAME);
    xvg_line = Malloc(char, MAX_LINE_FILE);
    aux_str = Malloc(char, 60);
    sprintf(file_name,"plot_nfront_%d_%d.xvg",*front, *ger);
    char *fname = path_join_file(local_execute,file_name);
    FILE *xvg_file = open_file(fname,fWRITE);
    fprintf (xvg_file,"#Obj1\tObj2\tValue_Ind_BEFORE_SORTED\n");
    for (int i=0; i < *size; i++){
        for (int obj = 0; obj < *num_obj; obj++){
            v_aux = get_displayed_value_of_objective(solutions_sorted, &i, &obj, 
                fitness_energies);            
            sprintf(aux_str,"%f",v_aux);
            if (obj == 0){
                strcpy(xvg_line, aux_str);
            }else{
                strcat(xvg_line, aux_str);
            }
            strcat(xvg_line, "\t");
        }
        //Obtaing individual when it was NOT sorted
        ind_id_before_sorted = get_nsga2_solution_index_by_objective_value(solutions_full, size_full,
            obj_get_index, &solutions_sorted[i].obj_values[*obj_get_index]) +1;
        sprintf(aux_str,"%d",ind_id_before_sorted);
        strcat(xvg_line, aux_str);
        fprintf (xvg_file,"%s\n", xvg_line);
    }
    fclose(xvg_file);
    free(aux_str);
    free(xvg_line);
    free(file_name);
    free(fname);
}


/** Saving files of initial population to analysis
*/
void saving_file_of_initial_population(const ea_nsga2_t *solutions_p, 
    const int *size, const input_parameters_t *in_para){

    solution_t *solutions_non_dominated = NULL;
    char *pop_non_dominated = NULL;
    int num_obj, number_of_non_dominated, front_non_dominated;    
    int ger;

    num_obj = solutions_p[0].sol->num_obj;

/**** Setting solutions to save information to NON-DOMINATED ***/
    ger = 0;
    front_non_dominated = 0;
    number_of_non_dominated = compute_how_many_front(solutions_p, size, &front_non_dominated);
    //Setting objectivies    
    solutions_non_dominated = allocate_solution(&number_of_non_dominated, &num_obj);
    set_nsga2_solution_in_solution(solutions_non_dominated, solutions_p, &number_of_non_dominated);
    //Setting proteins    
    copy_nsga2_solutions2solution(solutions_non_dominated, solutions_p, &number_of_non_dominated);    

    //Sorting solutions by first objective
    index_obj_compare_objective_solution = 0;
    qsort (solutions_non_dominated, number_of_non_dominated , sizeof (solution_t),
                      compare_objective_solution);

    //Saving data of non-dominated
    pop_non_dominated = Malloc(char, MAX_FILE_NAME);
    sprintf(pop_non_dominated,"pop_NON_DOMINATED_%d.pdb",ger);
    save_solution_population_file(solutions_non_dominated, in_para->path_local_execute, pop_non_dominated, 
        &number_of_non_dominated);
    build_objective_files_non_dominated(solutions_non_dominated, &ger, 
        &number_of_non_dominated, &num_obj, in_para);
    build_plot_xvg_file_initial_population(solutions_non_dominated, &number_of_non_dominated, &num_obj,
        in_para->path_local_execute, &front_non_dominated, &ger,
        in_para->fitness_energies, &index_obj_compare_objective_solution);
    //Desalocate protein in solution
    for (int s = 0; s < number_of_non_dominated; s++){
        protein_t* prot_aux_dest = (protein_t*) solutions_non_dominated[s].representation;
        if (prot_aux_dest != NULL){
            desAllocate_pdbatom(prot_aux_dest->p_atoms);
            prot_aux_dest->p_topol = NULL;
            desallocateProtein(prot_aux_dest, &number_of_non_dominated);
            solutions_non_dominated[s].representation = NULL;            
        }        
    }
    desallocate_solution(solutions_non_dominated, &number_of_non_dominated);
    free(pop_non_dominated);    
/**** FINISHED creating file to NON-DOMINATED **/

/**** Setting solutions to save information from solutions_p ***/   
    //Saving data of solution_p
    int one = 1;
    solution_t *solutions_p_AUX = NULL;    
    solutions_p_AUX = allocate_solution(&in_para->size_population, &num_obj);    
    //Setting objectivies    
    set_nsga2_solution_in_solution(solutions_p_AUX, solutions_p, &in_para->size_population);
    //Saving fitness values    
    build_fitness_files(solutions_p_AUX, &ger, &in_para->size_population);    
    desallocate_solution(solutions_p_AUX, &in_para->size_population);
/**** FINISHED creating file to solutions_p **/

}

void saving_file_to_generation_analysis(const ea_nsga2_t *solutions_rt, 
    const int *size_RT, const int *ger, const ea_nsga2_t *solutions_p,
    const input_parameters_t *in_para){

    solution_t *solutions_non_dominated = NULL; //*solutions,
    char *pop_non_dominated = NULL; //*pop_RT_file_name, *pop_p_file_name;
    int num_obj, number_of_non_dominated, front_non_dominated;    

    num_obj = solutions_rt[0].sol->num_obj;

/**** Setting solutions to save information to NON-DOMINATED ***/
    front_non_dominated = 0;
    number_of_non_dominated = compute_how_many_front(solutions_rt, size_RT, &front_non_dominated);
    //Setting objectivies    
    solutions_non_dominated = allocate_solution(&number_of_non_dominated, &num_obj);
    set_nsga2_solution_in_solution(solutions_non_dominated, solutions_rt, &number_of_non_dominated);
    //Setting proteins    
    copy_nsga2_solutions2solution(solutions_non_dominated, solutions_rt, &number_of_non_dominated);    

    //Sorting solutions by first objective
    index_obj_compare_objective_solution = 0;
    qsort (solutions_non_dominated, number_of_non_dominated , sizeof (solution_t),
                      compare_objective_solution);

    //Saving data of non-dominated
    pop_non_dominated = Malloc(char, MAX_FILE_NAME);
    sprintf(pop_non_dominated,"pop_NON_DOMINATED_%d.pdb",*ger);
    save_solution_population_file(solutions_non_dominated, in_para->path_local_execute, pop_non_dominated, 
        &number_of_non_dominated);
    build_objective_files_non_dominated(solutions_non_dominated, ger, 
        &number_of_non_dominated, &num_obj, in_para);
    build_plot_xvg_file(solutions_non_dominated, &number_of_non_dominated, &num_obj,
        in_para->path_local_execute, solutions_rt, size_RT, &front_non_dominated, ger,
        in_para->fitness_energies, &index_obj_compare_objective_solution);
    //Desalocate protein in solution
    for (int s = 0; s < number_of_non_dominated; s++){
        protein_t* prot_aux_dest = (protein_t*) solutions_non_dominated[s].representation;
        if (prot_aux_dest != NULL){
            desAllocate_pdbatom(prot_aux_dest->p_atoms);
            prot_aux_dest->p_topol = NULL;
            desallocateProtein(prot_aux_dest, &number_of_non_dominated);
            solutions_non_dominated[s].representation = NULL;            
        }        
    }
    desallocate_solution(solutions_non_dominated, &number_of_non_dominated);
    free(pop_non_dominated);    
/**** FINISHED creating file to NON-DOMINATED **/

/**** Setting solutions to save information from solutions_p ***/   
    //Saving data of solution_p
    int one = 1;
    solution_t *solutions_p_AUX = NULL;
    char *pop_p_file_name = NULL;
    pop_p_file_name = Malloc(char, MAX_FILE_NAME);
    solutions_p_AUX = allocate_solution(&in_para->size_population, &num_obj);    
    //Setting objectivies    
    set_nsga2_solution_in_solution(solutions_p_AUX, solutions_p, &in_para->size_population);
    //Setting proteins    
    copy_nsga2_solutions2solution(solutions_p_AUX, solutions_p, &in_para->size_population);    
    //Saving fitness values    
    build_fitness_files(solutions_p_AUX, ger, &in_para->size_population);    
    sprintf(pop_p_file_name,"pop_%d.pdb",*ger);
    save_solution_population_file(solutions_p_AUX, in_para->path_local_execute, pop_p_file_name, 
        &in_para->size_population);    
    //Desalocate protein in solution
    for (int s = 0; s < in_para->size_population; s++){
        protein_t* prot_aux_dest = (protein_t*) solutions_p_AUX[s].representation;
        if (prot_aux_dest != NULL){
            desAllocate_pdbatom(prot_aux_dest->p_atoms);
            prot_aux_dest->p_topol = NULL;
            desallocateProtein(prot_aux_dest, &one);
            solutions_p_AUX[s].representation = NULL;            
        }        
    }    
    desallocate_solution(solutions_p_AUX, &in_para->size_population);
    free(pop_p_file_name);
/**** FINISHED creating file to solutions_p **/


/**** Setting solutions to save information from solutions_rt ***/
/*    
    //Setting objectivies
    solutions = allocate_solution(size_RT, &num_obj);    
    set_nsga2_solution_in_solution(solutions, solutions_rt, size_RT);

    //Saving data of solutions_rt
    pop_RT_file_name = Malloc(char, MAX_FILE_NAME);
    sprintf(pop_RT_file_name,"pop_RT_%d.pdb",*ger);
    save_nsga2_population_file(solutions_rt, in_para->path_local_execute, pop_RT_file_name, 
        size_RT);
    build_fitness_files(solutions, ger, size_RT);
    free(pop_non_dominated);
    free(pop_RT_file_name);
    desallocate_solution(solutions, size_RT);    
*/    
/**** FINISHED creating file to solutions_rt **/
}


/** Copies survivors to next generation. 
* solutions_p stores the individuals which will survive to next generation
* solutions_rt contains all solutions of current generation
* size is the size to nsga2_solutions_p
* Important: nsga2_solutions_rt is sorted by front and crowding distance.
* Therefore, the first N individuals are the best solutions.
*/
void copy_survivors_to_next_generation(ea_nsga2_t *solutions_p, 
            const ea_nsga2_t *solutions_rt, const int *size){
    for(int i = 0; i < *size; i++){
        copy_nsga2_solution(&solutions_p[i], &solutions_rt[i]);
    }
}

_2PG_CARTESIAN_EXPORT
int ea_nsga2(const input_parameters_t *in_para){
    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    
    protein_t *population_p; // main population
    /* pop_new is genetared by reproduce function. 
     * It will  set in solution_q
    */
    protein_t *pop_new;

    ea_nsga2_t *nsga2_solutions_p; // main solution
    ea_nsga2_t *nsga2_solutions_q; // new solutions that stores child of nsga2_solutions_p
    ea_nsga2_t *nsga2_solutions_rt; // solution that join p and q

    int generation;
    int size_RT;

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

    //Allocating PDB ATOMS of population_p
    population_p = allocateProtein(&in_para->size_population);    

    //Loading initial population
    load_initial_population_file(population_p, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);

    //Building Global topology since all models have same topology
    int one = 1;
    protein_ref_global = allocateProtein(&one);
    protein_ref_global->p_atoms = allocate_pdbatom(&population_p[0].p_topol->numatom);
    copy_protein_atoms(protein_ref_global, &population_p[0]);
    protein_ref_global->p_topol = allocateTop_Global(&population_p[0].p_topol->numres, 
            &population_p[0].p_topol->numatom);    
    build_topology_individual(protein_ref_global);
    nsga2_global_topol = protein_ref_global->p_topol;

    //Setting population_new
    pop_new = allocateProtein(&in_para->size_population);
    //Allocating protein at pop_new because copy_portein does not allocated anymore    
    for (int i = 0; i < in_para->size_population; i++){
        pop_new[i].p_atoms = allocate_pdbatom(&nsga2_global_topol->numatom);
        pop_new[i].p_topol = population_p[0].p_topol;
    }       
    copy_protein_population(pop_new, population_p, &in_para->size_population);
    initialize_protein_population_atoms(pop_new, &in_para->size_population);

/*
    //Saving topology of population 
    char *prefix = Malloc(char,10);
    strcpy(prefix, "prot");
    save_topology_population(population_p, &in_para->size_population, 
    in_para->path_local_execute, prefix);    
    free(prefix);
*/    
/**************** STARTING NSGA-II Algorithm *************************/
    display_msg("Starting NSGA-II Algorithm \n");
    initialize_algorithm_execution(primary_sequence, in_para);
    init_gromacs_execution();

    //Setting solutions nsga2_solutions_p   
    nsga2_solutions_p = allocate_nsga2(in_para);
    //Setting solutions nsga2_solutions_q
    nsga2_solutions_q = allocate_nsga2(in_para);
    //Setting solutions nsga2_solutions_rt
    size_RT = 2*in_para->size_population;  
    nsga2_solutions_rt = allocate_nsga2_RT(&size_RT, in_para);

    //Setting reference of proteins to solution 
    set_proteins2nsga2_solutions(nsga2_solutions_p, population_p, 
    	&in_para->size_population);
    
    //Computing objectives of solutions p with GROMACS
    nsga2_compute_objectives_with_gromacs(nsga2_solutions_p, 
    	&in_para->size_population, in_para);

    //Setting fronts to nsga2_solutions_p
    set_front_solution_p(nsga2_solutions_p, &in_para->size_population);

    //Saving information of initiak population to analysis
    saving_file_of_initial_population(nsga2_solutions_p, &in_para->size_population, in_para);

    // Main Looop of NSGA-II Algorithm 
    int started_generation = get_started_generation(&in_para->started_generation);
    for (int g = started_generation; g < (in_para->number_generation+started_generation); g++){
        if (started_generation == 0){
            generation = g +1; //used to messages. It is started to 1.        
        }else{
            generation = g; //since g is greater than zero
        }        
        
        //Reproduce nsga2_solutions_p to generate pop_new
        reproduce_protein(pop_new, nsga2_solutions_p, in_para);

        //Setting reference of new proteins to solution_q 
        set_proteins2nsga2_solutions(nsga2_solutions_q, pop_new, 
            &in_para->size_population);

        //Computing objectives of solutions q with GROMACS
        nsga2_compute_objectives_with_gromacs(nsga2_solutions_q, 
            &in_para->size_population, in_para);

        //Set solution_p and solution_q in solution_rt
        set_soluton_rt(nsga2_solutions_rt, 
            nsga2_solutions_p, nsga2_solutions_q, 
            &in_para->size_population);

        //Set and Sorting population rt based on Dominance and Crowding Distance
        set_dominance_and_crowding_distance_in_soluton_rt(nsga2_solutions_rt, 
            &size_RT, in_para);

        //Coping to next generation
        copy_survivors_to_next_generation(nsga2_solutions_p, 
            nsga2_solutions_rt, &in_para->size_population);

        //Creating files for generation analysis
        saving_file_to_generation_analysis(nsga2_solutions_rt, &size_RT,
            &generation, nsga2_solutions_p, in_para);
    }
	finish_gromacs_execution();
/**************** FINISHED NSGA-II Algorithm *************************/

    desAllocateTop_Global(nsga2_global_topol);
    desAllocate_pdbatom(protein_ref_global->p_atoms);
    desallocateProtein(protein_ref_global, &one);    
    desallocate_solution_nsga2_RT(nsga2_solutions_rt, &size_RT);
    desallocate_solution_nsga2(nsga2_solutions_p, &in_para->size_population);
    desallocate_solution_nsga2(nsga2_solutions_q, &in_para->size_population);    
    desallocateProtein(pop_new, &in_para->size_population);        
    desallocateProtein(population_p, &in_para->size_population);
    desallocate_primary_seq(primary_sequence);
 
    return 0;    
	
}
