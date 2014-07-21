#include <stdlib.h>
#include <string.h>

#ifndef WIN32
#include <unistd.h>
#include <dirent.h>
#endif

#include <sys/stat.h>
#include <sys/types.h>

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
#include "owner_file_analysis.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

#define MAX_LINE_SOLUTION_FILE 300

_2PG_CARTESIAN_EXPORT
owner_file_t* allocate_file_t(const int *num_files, const int *num_obj){
	owner_file_t* file_names_aux = NULL;
	file_names_aux = Malloc(owner_file_t, *num_files);
	for (int i = 0; i < *num_files; i++){
		file_names_aux[i].file_name = Malloc(char, MAX_FILE_NAME);
		file_names_aux[i].obj_values = Malloc(double, *num_obj);
	}		
	return file_names_aux;
}

_2PG_CARTESIAN_EXPORT
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

int compare_first_value(const void *x, const void *y){
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

int compare_front(const void *x, const void *y){
    int fx, fy;
    fx = ((owner_file_t *)x)->front;
    fy = ((owner_file_t *)y)->front;
    if (fx > fy){
        return 1;
    }else {
        return 0;
    }
}

int compare_dominated(const void *x, const void *y){
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


void save_analysis_files_no_objectives(const owner_file_t *solutions_f, const int *size, const int *numobj, const char *column_1, const char *column_2 ){
	
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
    strcpy(c_obj, column_1);
    strcpy(c_obj_2, column_2);

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
		    		solutions_f[s].obj_values[ob1],
		    		solutions_f[s].obj_values[ob2],
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
		    		temp_aux[j].obj_values[ob1],
		    		temp_aux[j].obj_values[ob2],
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



owner_file_t * loading_owner_file_solution(const int *num_solutions_r, 	const int *numobj_r, const char *path_file_name){
	FILE * solution_file;	
	char *line;
	char *line_splited;
	int sol;
	owner_file_t * solutions_aux;

	line = Malloc(char, MAX_LINE_SOLUTION_FILE);

	//Alocating Solution
	solutions_aux = allocate_file_t(num_solutions_r, numobj_r);

	//Reading file and set values of objective
	sol = -1;
	solution_file = open_file(path_file_name, fREAD);
	//Removing first line that is collumn
	fgets(line,MAX_LINE_SOLUTION_FILE,solution_file);
	while ( fgets(line,MAX_LINE_SOLUTION_FILE,solution_file) != NULL){
		sol = sol + 1;
		//Obtaing index collumn
		line_splited = strtok (line,"\t");
		strcpy(solutions_aux[sol].file_name, line_splited);
		//Setting number of objectives		
		for (int ob = 0; ob < *numobj_r; ob++){
			line_splited = strtok (NULL,"\t");
			trim(line_splited);
			solutions_aux[sol].obj_values[ob] = str2double(line_splited);
		}					
	}	
	fclose(solution_file);

	free(line);

	return solutions_aux;
}
