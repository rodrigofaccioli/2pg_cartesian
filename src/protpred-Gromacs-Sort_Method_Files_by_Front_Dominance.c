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


typedef struct sowner_file{	
	char *file_name;
	double *obj_values;
	int front;
}owner_file_t;

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


void save_front(const owner_file_t *solutions_f, const int *size, const int *numobj, const type_fitness_energies_t *fitness_energies){
	FILE *d_file=NULL;
	char *file_name=NULL;
	char *line=NULL;
	char *aux_str=NULL;

	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "front_methods_output.txt");	
	line = Malloc(char, MAX_LINE_FILE);
    aux_str = Malloc(char, 60);

	d_file = open_file(file_name, fWRITE);
	fprintf(d_file, "#Method\tValue1\tValue2\tFront\n");
	for (int s = 0; s < *size; s++){
		strcpy(line, solutions_f[s].file_name);
		strcat(line,"\t");
		for (int ob = 0; ob < *numobj; ob++){
			sprintf(aux_str,"%f", get_displayed_value_of_objective_from_owner_file_t(solutions_f, &s, &ob, 
                fitness_energies));
			strcat(line,aux_str);
			strcat(line,"\t");
		}
		int2str(aux_str, &solutions_f[s].front);
		strcat(line,aux_str);
		strcat(line,"\t");
		fprintf(d_file, "%s\n",line);
	}
	fclose(d_file);
	free(aux_str);
	free(line);
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
	display_msg("Run of pdb2gmx for pattern of atom names  \n");
	for (int ind = 0; ind < num_files; ind ++){
		call_pdb2gmx_for_pattern_atom_names(file_names[ind].file_name, in_param.path_local_execute, in_param.path_gromacs_programs, in_param.force_field);
		clean_gromacs_simulation(in_param.path_local_execute);
	}	
	// FINISHED RUN PDB2GMX 

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
    //Sorting by front
    qsort(file_names, num_files,  sizeof (owner_file_t), compare_front);

    //Saving file
	save_front(file_names, &num_files, &in_param.number_fitness, in_param.fitness_energies);
/**************** FINISHED FINAL RESULTS *************************/

	desalocate_file_t(file_names, &num_files);
    desallocate_solution(solutions_p, &num_files);
    desallocateProtein(population_p, &num_files);	
	deAllocateload_parameters(&in_param);
	display_msg("Done !!! \n");
	return 0;
}