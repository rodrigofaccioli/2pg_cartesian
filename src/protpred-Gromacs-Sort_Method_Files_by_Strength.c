#ifdef _WIN32
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent_2pgwin.h>
#include <time.h>
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
//#include "ea_spea2.h"
#include "dominancep.h"
#include "spea2.h"
#include "owner_file_analysis.h"
#endif

#ifdef linux
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <dirent.h>
#include <time.h>
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
//#include "ea_spea2.h"
#include "dominancep.h"
#include "spea2.h"
#include "owner_file_analysis.h"
#endif

int main(int argc, char *argv[]){
	input_parameters_t in_param;
	display_msg("Sort by spea2 - Reading the configure file \n");
	load_parameters_from_file(&in_param,argv[1]);

	/*Represents extension of files that will be looked in directory*/
	char *ext = NULL;
	/*Represents how many files there is in directory*/
	int num_files;
	/*Represents file names that are in directory and have target extension*/
	owner_file_t *file_names = NULL;

	protein_t *population_p= NULL; // main population of protein
	solution_t *solution_p= NULL; // main solution
	dominancep_t *dominancep = NULL; // main of dominancep
	spea2_t *spea2 = NULL; // main of spea2
	


/**************** START GETTING FILE NAMES *************************/
	ext = Malloc(char, 4);
	strcpy(ext, "pdb");
	num_files = how_many_files_directory_by_extension(in_param.path_local_execute, ext);
	file_names = allocate_file_t(&num_files, &in_param.number_fitness);
	insert_files_directory_by_extension(file_names, in_param.path_local_execute, ext);
	//free(ext);
/**************** FINISHED GETTING FILE NAMES *************************/


/**************** START GETTING THE OBJECTIVES *************************/	
	init_gromacs_execution();
/*
	// RUN PDB2GMX IN ALL FILES FOR PATTERN OF ATOM NAMES 
	display_msg("Run of pdb2gmx for pattern of atom names\n");
	for (int ind = 0; ind < num_files; ind ++){
		call_pdb2gmx_for_pattern_atom_names(file_names[ind].file_name, in_param.path_local_execute, in_param.path_gromacs_programs, in_param.force_field);
		clean_gromacs_simulation(in_param.path_local_execute);
	}	
	// FINISHED RUN PDB2GMX 
*/
	population_p = allocateProtein(&num_files);
	solution_p = allocate_solution(&num_files, &in_param.number_fitness);        
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
		solution_p[ind].representation = &population_p[ind];
		get_gromacs_objectives_of_solution(&solution_p[ind], &in_param, &ind);
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
	//Setting identification
	for (int i = 0; i < num_files; i++){
	   solution_p[i].ID = i+1;
	}

    //Setting dominancep and spea2
	dominancep = allocate_dominancep(&num_files);
	spea2 = allocate_spea2(&num_files);
//	printf("set_dominancep\n");
/*
	solution_p[0].obj_values[0] = 5;
	solution_p[0].obj_values[1] = 3;
	solution_p[1].obj_values[0] = 2;
	solution_p[1].obj_values[1] = 8;
	solution_p[2].obj_values[0] = 3;
	solution_p[2].obj_values[1] = 2;
	solution_p[3].obj_values[0] = 5;
	solution_p[3].obj_values[1] = 1;
	solution_p[4].obj_values[0] = 8;
	solution_p[4].obj_values[1] = 7;
	solution_p[5].obj_values[0] = 1;
	solution_p[5].obj_values[1] = 2;
*/
//	clock_t tic = clock();	
	set_dominancep(dominancep, solution_p, &num_files);
//	printf("show_dominance\n");
	//show_dominancep(dominancep, &num_files);

	initialize_spea2(spea2, dominancep, &num_files, file_names);
	//printf("set_raw_fitness from spea2\n\n");
	set_raw_fitness(spea2, dominancep, &num_files);
	//printf("set_density from spea2\n\n");
	set_density(spea2, &num_files, &num_files);
	//printf("set_fitness from spea2\n\n");
	set_fitness(spea2, &num_files);
	//printf("order from spea2\n\n");
	order(spea2, &num_files);
//	clock_t toc = clock();
//	printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);	        
	save_spea2_v1(spea2, &num_files);
        save_spea2_v2(spea2, &num_files);

	

	desallocate_dominancep(dominancep, &num_files);
	desallocate_spea2(spea2);
/**************** FINISHED GETTING FRONT *************************/

/**************** START GETTING FINAL RESULTS *************************/
	/*
    //Sorting solutions
    sorting_solutions_by_front_dominance(file_names, &num_files, &in_param.number_fitness);

    //Saving file
	save_analysis_files(file_names, &num_files, &in_param.number_fitness, in_param.fitness_energies);
*/
/**************** FINISHED FINAL RESULTS *************************/

	desalocate_file_t(file_names, &num_files);
        desallocate_solution(solution_p, &num_files);
        desallocateProtein(population_p, &num_files);	
	deAllocateload_parameters(&in_param);
	display_msg("Done !!! \n");
	return 0;
}
