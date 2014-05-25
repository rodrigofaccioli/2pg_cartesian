#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>

#include "defines.h"
#include "load_parameters.h"
#include "messages.h"
#include "defines.h"
#include "enums.h"
#include "ea_mono.h"
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
#include "ea_nsga2.h"
#include "dominance.h"
#include "analysis_types.h"


owner_file_t* allocate_file_t(const int *num_files, const int *num_obj){
	owner_file_t* file_names_aux = NULL;
	file_names_aux = Malloc(owner_file_t, *num_files);
	for (int i = 0; i < *num_files; i++){
		file_names_aux[i].file_name = Malloc(char, MAX_FILE_NAME);
		file_names_aux[i].obj_values = Malloc(double, *num_obj);
	}		
	return file_names_aux;
}

void desalocate_file_t(owner_file_t*file_names, const int *num_files){
	for (int i = 0; i < *num_files; i++){
		free(file_names[i].file_name);
		free(file_names[i].obj_values);
	}	
	free(file_names);
	file_names = NULL;	
}

void copy_file_owner(owner_file_t* dest, const owner_file_t*source, const int *num){
	strcpy(dest->file_name, source->file_name);
	for (int i = 0; i < *num; i++){
		dest->obj_values[i] = source->obj_values[i];
	}
	dest->front = source->front;
	dest->number_solutions_are_dominated = source->number_solutions_are_dominated;
	dest->ranking = source->ranking;

}

const char *get_filename_ext(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if(!dot || dot == filename) return "";
    return dot + 1;
}

int how_many_files_directory_by_extension(const char *path, const char *ext){
	const char *ext_aux = NULL;	
	int r = 0;
	DIR *dir = opendir(path);
	struct dirent *ent;
	while (ent = readdir(dir)) {
		ext_aux = get_filename_ext(ent->d_name);
		if (strcmp (ext_aux, ext) == 0){
			r++;
		}	
	}
	closedir(dir);
	return r;
}

void insert_files_directory_by_extension(owner_file_t *file_names, const char *path, const char *ext){
	const char *ext_aux = NULL;
	int index = -1;

	DIR *dir = opendir(path);
	struct dirent *ent;
	while (ent = readdir(dir)) {
		ext_aux = get_filename_ext(ent->d_name);
		if (strcmp (ext_aux, ext) == 0){
			index = index +1;
			strcpy(file_names[index].file_name, ent->d_name);
		}	
	}
	closedir(dir);
}

/** Returns the value of an specific objective */
double get_specific_objective_from_owner_file_t(const owner_file_t *sol, const int *obj){
	return sol->obj_values[*obj];
}


/** Returns the oposite value of an specific objective */
double get_oposite_specific_objective_from_owner_file_t(const owner_file_t *sol, const int *obj){
	/* objectives that must be maximized when they are obtained these values
	 * are multiplied by -1 because 2PG works considerating minimization.
	 * However, when they will be stored, they must be written in original
	 * value 
	 */	
	return sol->obj_values[*obj] * (-1);
}


/** Returns the displyed value of an specific objective */
double get_displayed_value_of_objective_from_owner_file_t(const owner_file_t *sol, 
	const int *index, const int *obj, 
	const type_fitness_energies_t *fitness_energies){
	if ( (fitness_energies[*obj] == fit_hbond) ||
		(fitness_energies[*obj] == fit_hydrophilic)||
		(fitness_energies[*obj] == fit_hbond_main)	||
		(fitness_energies[*obj] == fit_stride_total) ||
		(fitness_energies[*obj] == fit_stride_helix)	||
		(fitness_energies[*obj] == fit_stride_beta) ) {
		return get_oposite_specific_objective_from_owner_file_t(&sol[*index],obj);
	}else{
		return get_specific_objective_from_owner_file_t(&sol[*index],obj);
	}
}

/** Computes how many individuals there is in front 
*/
int compute_how_many_front_file_t(const owner_file_t * solutions, 
            const int *size, const int *front_ref){
    int num = 0;
    for (int i = 0; i < *size; i++){
        if (solutions[i].front == *front_ref){
            num = num + 1;
        }
    }
    return num;
}

static int compare_first_value(const void *x, const void *y){
    int fx, fy;
    fx = ((owner_file_t *)x)->obj_values[0];
    fy = ((owner_file_t *)y)->obj_values[0];
    if (fx > fy){
        return 1;
    }else {
        return 0;
    }
}

void save_analysis_files(const owner_file_t *solutions_f, const int *size, const int *numobj, const type_fitness_energies_t *fitness_energies){
	
	FILE *d_file=NULL;
	char *file_name=NULL;	
	char *line_f=NULL;
	char *aux_str=NULL;
	char *c_front=NULL;	
	char *c_obj=NULL;	
	char *c_obj_2=NULL;
	owner_file_t *temp_aux=NULL;
	int ob1, ob2 = -1;
	file_name = Malloc(char, 100);	
	line_f = Malloc(char, MAX_LINE_FILE);
    aux_str = Malloc(char, 60);
    c_front = Malloc(char, 10);
    c_obj = Malloc(char, MAX_RANDOM_STRING);
    c_obj_2 = Malloc(char, MAX_RANDOM_STRING);
    ob1 = 0;
    ob2 = 1;
    int front;
    int i;
    int how_many_front;
    int ind_temp;

    //Getting name of objectives
    type_fitness_energies2str(c_obj, &fitness_energies[0]);
    type_fitness_energies2str(c_obj_2, &fitness_energies[1]);

 /**************** STARTING FRONT FILES *************************************/
	front = 0;    
    strcpy(file_name, "all_front_");
    strcat(file_name,c_obj);
    strcat(file_name,"_");    	    
    strcat(file_name,c_obj_2);
    strcat(file_name,".front");    
    d_file = open_file(file_name, fWRITE);
	fprintf(d_file, "#Ranking represents classification based on all solutions\n");    
	fprintf(d_file, "#Front represents a classification based on dominance critera\n");
	fprintf(d_file, "#Dominated represents number of solutions are dominated by solution\n");					
	sprintf(line_f, "#%s is the first value used to apply dominance critera\n", c_obj);			
	fprintf(d_file, line_f);
	sprintf(line_f, "#%s is the second value used to apply dominance critera\n", c_obj_2);
	fprintf(d_file, line_f);
	fprintf(d_file, "#Method represents name of method which is a solution\n");	
	fprintf(d_file, "#\n");
	fprintf(d_file, "#\n");
	sprintf(line_f, "#Ranking  Front\tDominated\t%s\t\t%s\t\tMethod\n", c_obj, c_obj_2 );
	fprintf(d_file, line_f);
    for (int s = 0; s < *size; s++){
	   	fprintf(d_file, "%5d\t  %3d\t%5d\t\t%.10e\t%.10e\t%-40.40s\n", 
		    		solutions_f[s].ranking,
		    		solutions_f[s].front,
		    		solutions_f[s].number_solutions_are_dominated,
					get_displayed_value_of_objective_from_owner_file_t(solutions_f, &s, &ob1, fitness_energies), 
					get_displayed_value_of_objective_from_owner_file_t(solutions_f, &s, &ob2, fitness_energies), 
					solutions_f[s].file_name);
    }
    fclose(d_file);
/**************** FINISHED FRONT FILES *************************************/

/**************** STARTING XVG FILES *************************************/
    front = 0;
    i = 0;
    while (i < *size){
    	    strcpy(file_name, "plot_front_");
    	    int2str(c_front, &front);
    	    strcat(file_name,c_front);    	    
    	    strcat(file_name,"_");
    	    strcat(file_name,c_obj);
    	    strcat(file_name,"_");    	    
    	    strcat(file_name,c_obj_2);
		    strcat(file_name,".xvg");
		    d_file = open_file(file_name, fWRITE);
			sprintf(line_f, "#%s is the first value used to apply dominance critera\n", c_obj);			
			fprintf(d_file, line_f);
			sprintf(line_f, "#%s is the second value used to apply dominance critera\n", c_obj_2);
			fprintf(d_file, line_f);		    
			fprintf(d_file, "#Front represents a classification based on dominance critera\n");
			fprintf(d_file, "#Dominated represents number of solutions are dominated by solution\n");			
			fprintf(d_file, "#Ranking represents classification based on all solutions\n");
			fprintf(d_file, "#\n");
			fprintf(d_file, "#\n");
			sprintf(line_f, "#%s\t\t%s\t\tFront\tDominated\tRanking\n", c_obj, c_obj_2 );
			fprintf(d_file, line_f);
			//Preparing temporary solutions
			how_many_front = compute_how_many_front_file_t(solutions_f, size, &front);
			temp_aux = allocate_file_t(&how_many_front, numobj);
			ind_temp = 0;
		    while (solutions_f[i].front == front){
		    	//Coping solutions to temporary that will be sorted
		    	copy_file_owner(&temp_aux[ind_temp], &solutions_f[i], numobj);
		    	ind_temp = ind_temp + 1;
		    	i = i + 1;		    	
		    }
		    //Sorting temporary solutions by first value
		    qsort(temp_aux, how_many_front, sizeof(owner_file_t), compare_first_value);
		    for (int j = 0; j < how_many_front; j++){
		    	fprintf(d_file, "%.10e \t%.10e\t%3d\t%5d\t%10d\n", 
					get_displayed_value_of_objective_from_owner_file_t(temp_aux, &j, &ob1, fitness_energies), 
					get_displayed_value_of_objective_from_owner_file_t(temp_aux, &j, &ob2, fitness_energies), 
					temp_aux[j].front,
					temp_aux[j].number_solutions_are_dominated,
					temp_aux[j].ranking);
		    }
		    fclose(d_file); //close file of front		    
		    desalocate_file_t(temp_aux, &how_many_front); //desalocate temporary
		    front = front + 1; //Next front
    }
/**************** FINISHED XVG FILES *************************************/    
    free(c_obj_2);
    free(c_obj);
	free(c_front);
	free(line_f);
	free(aux_str);
	free(file_name);		
}

static int compare_front(const void *x, const void *y){
    int fx, fy;
    fx = ((owner_file_t *)x)->front;
    fy = ((owner_file_t *)y)->front;
    if (fx > fy){
        return 1;
    }else {
        return 0;
    }
}

static int compare_dominated(const void *x, const void *y){
    int fx, fy;
    fx = ((owner_file_t *)x)->number_solutions_are_dominated;
    fy = ((owner_file_t *)y)->number_solutions_are_dominated;
    if (fx < fy){
        return 1;
    }else {
        return 0;
    }
}

/** In this function is assigned a value based on a sorting
*/
void set_global_ranking(owner_file_t *file_names, const int *size){

	for (int i = 0; i < *size; i++){
		file_names[i].ranking = i+1;
	}

}

/** Solutions are sorted based on front in crescent order . 
* After, each front is sorted by number of dominated solutions 
*/
void sorting_solutions_by_front_dominance(owner_file_t *file_names, const int *size, const int *numobj){	
	int how_many_front;	
	int front;
	int start_index;
	int total_ind;
	int j;
	owner_file_t *aux=NULL;		

	//Sorting by Front
	qsort(file_names, *size,  sizeof (owner_file_t), compare_front);

	//Sorting by dominated in each front
	front = 0;	
	start_index = 0;
	how_many_front = compute_how_many_front_file_t(file_names, size, &front);
	while (how_many_front > 0){
		aux = allocate_file_t(&how_many_front, numobj);
		total_ind = start_index + how_many_front;
		j = 0;
		for (int index_f = start_index; index_f < total_ind; index_f++){
			copy_file_owner(&aux[j], &file_names[index_f], numobj);
			j = j + 1;
		}
		qsort(aux, how_many_front,  sizeof (owner_file_t), compare_dominated);		
		j = 0;
		for (int index_f = start_index; index_f < total_ind; index_f++){
			copy_file_owner(&file_names[index_f], &aux[j], numobj);
			j = j + 1;
		}		
		desalocate_file_t(aux, &how_many_front);
		start_index = start_index + how_many_front;
		front = front + 1;		
		how_many_front = compute_how_many_front_file_t(file_names, size, &front);
	}
	set_global_ranking(file_names, size);
}

int main(int argc, char *argv[]){
	input_parameters_t in_param;
	display_msg("Reading the configure file \n");
	load_parameters_from_file(&in_param,argv[1]);

	/*Represents extension of files that will be looked in directory*/
	char *ext = NULL;
	/*Represents how many files there is in directory*/
	int num_files;
	/*Represents file names that are in directory and have target extension*/
	owner_file_t *file_names = NULL;

	protein_t *population_p= NULL; // main population of protein
	solution_t *solutions_p= NULL; // main solution
	dominance_t *dominance= NULL; // main of dominance
	ea_nsga2_t *nsga2_solutions_p= NULL; // main of front


/**************** START GETTING FILE NAMES *************************/
	ext = Malloc(char, 4);
	strcpy(ext, "pdb");
	num_files = how_many_files_directory_by_extension(in_param.path_local_execute, ext);
	file_names = allocate_file_t(&num_files, &in_param.number_fitness);
	insert_files_directory_by_extension(file_names, in_param.path_local_execute, ext);
	free(ext);
/**************** FINISHED GETTING FILE NAMES *************************/


/**************** START GETTING THE OBJECTIVES *************************/	
	init_gromacs_execution();

	// RUN PDB2GMX IN ALL FILES FOR PATTERN OF ATOM NAMES 
	display_msg("Run of pdb2gmx for pattern of atom names\n");
	for (int ind = 0; ind < num_files; ind ++){
		call_pdb2gmx_for_pattern_atom_names(file_names[ind].file_name, in_param.path_local_execute, in_param.path_gromacs_programs, in_param.force_field);
		clean_gromacs_simulation(in_param.path_local_execute);
	}	
	// FINISHED RUN PDB2GMX 
/*
	// RUN MINIMIZATION OF ALL SOLUTIONS	
	display_msg("Run of minimization all solutions\n");	
	for (int ind = 0; ind < num_files; ind ++){
	   	build_tpr_file(file_names[ind].file_name, in_param.path_local_execute, in_param.path_gromacs_programs, 
	        	in_param.force_field, in_param.mdp_file_min);
		call_mdrun2minimization(file_names[ind].file_name, in_param.path_local_execute, in_param.path_gromacs_programs);
		clean_gromacs_simulation(in_param.path_local_execute);
	}
	// FINISHED MINIMIZATION OF ALL SOLUTIONS
*/
	population_p = allocateProtein(&num_files);
	solutions_p = allocate_solution(&num_files, &in_param.number_fitness);        
	for (int ind = 0; ind < num_files; ind ++){
		printf("Computing Objectives of %s\n", file_names[ind].file_name);
		char *path_PDB_file_name = path_join_file(in_param.path_local_execute, file_names[ind].file_name);
		int num_atom = get_num_atom(path_PDB_file_name);		
		population_p[ind].p_atoms = allocate_pdbatom(&num_atom);		
		load_pdb_file_without_num_atom(population_p[ind].p_atoms, NULL,	path_PDB_file_name);
		int num_res = get_number_residues_from_atom(population_p[ind].p_atoms, &num_atom);
		population_p[ind].p_topol = allocateTop_Global(&num_res, &num_atom);		
		renumerate_residue_number(population_p[ind].p_atoms, &num_atom);		
		build_topology_individual(&population_p[ind]);		
		solutions_p[ind].representation = &population_p[ind];
		get_gromacs_objectives_of_solution(&solutions_p[ind], &in_param, &ind);
		desAllocateTop_Global(population_p[ind].p_topol);
		desAllocate_pdbatom(population_p[ind].p_atoms);
		free(path_PDB_file_name);
	}	
	finish_gromacs_execution();
	//Deleting temporary files
	system("rm PROT_IND_*.pdb");
/**************** FINISHED GETTING THE OBJECTIVES *************************/

/**************** START GETTING FRONT *************************/
	in_param.size_population = num_files;
	nsga2_solutions_p = allocate_nsga2_without_allocation_of_representation(&in_param);	
	//Setting identification
	for (int i = 0; i < num_files; i++){
		nsga2_solutions_p[i].sol->ID = i+1;
	}

    //Setting dominance
	dominance = allocate_dominance(&num_files);
	set_dominance(dominance, solutions_p, &num_files);

	//Coping values of dominance
	for (int ind = 0; ind < num_files; ind++){
		// Indicates number of solutions that are dominated by me
		file_names[ind].number_solutions_are_dominated = dominance[ind].max_dominated;
	}	

	//Coping values of objective
	for (int ind = 0; ind < num_files; ind++){
		for (int ob = 0; ob < in_param.number_fitness; ob++)	
			nsga2_solutions_p[ind].sol->obj_values[ob] = solutions_p[ind].obj_values[ob];
	}

    //Setting front based on dominance concept
    compute_fronts(nsga2_solutions_p, dominance, &num_files);

    //Coping values from nsga2_solutions_p to owner_file_t
	for (int ind = 0; ind < num_files; ind++){
		//Coping objectives
		for (int ob = 0; ob < in_param.number_fitness; ob++){	
			file_names[ind].obj_values[ob] = nsga2_solutions_p[ind].sol->obj_values[ob];
		}
		file_names[ind].front = nsga2_solutions_p[ind].front;
	}
	desallocate_dominance(dominance, &num_files);
	desallocate_solution_nsga2(nsga2_solutions_p, &num_files);
/**************** FINISHED GETTING FRONT *************************/

/**************** START GETTING FINAL RESULTS *************************/
    //Sorting solutions
    sorting_solutions_by_front_dominance(file_names, &num_files, &in_param.number_fitness);

    //Saving file
	save_analysis_files(file_names, &num_files, &in_param.number_fitness, in_param.fitness_energies);
/**************** FINISHED FINAL RESULTS *************************/

	desalocate_file_t(file_names, &num_files);
    desallocate_solution(solutions_p, &num_files);
    desallocateProtein(population_p, &num_files);	
	deAllocateload_parameters(&in_param);
	display_msg("Done !!! \n");
	return 0;
}