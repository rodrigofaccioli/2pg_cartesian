#include <stdlib.h>
#include <string.h>

#include "dominance.h"
#include "load_parameters.h"
#include "messages.h"
#include "defines.h"
#include "solution.h"
#include "futil.h"
#include "string_owner.h"
#include "solutionio.h"
#include "ea_nsga2.h"
#include "owner_file_analysis.h"

void set_columns_of_solution_file(char *column_1, char *column_2, const char *path_file_name){
	char *line_splited;
	int num_obj = 0;
	FILE * solution_file;
	char *line;

	line = Malloc(char, 300);
	solution_file = open_file(path_file_name, fREAD);

	//getting first line to obtain number of objectivies
	fgets(line,300,solution_file);
 	line_splited = strtok (line,"\t");
  	while (line_splited != NULL){
    	num_obj = num_obj + 1;
    	if (num_obj == 2){
    		strcpy(column_1, line_splited);
    	}else if (num_obj == 3){
    		strcpy(column_2, line_splited);
    	}
    	line_splited = strtok(NULL, "\t");
  	}	
  	fclose(solution_file);
  	trim(column_1);
  	trim(column_2);
  	remove_character(column_1, '#');
  	remove_character(column_2, '\n');
  	free(line);

}


int main(int argc, char *argv[]){
	input_parameters_t in_param;
	//display_msg("Reading the configure file \n");
	//load_parameters_from_file(&in_param,argv[1]);

	dominance_t *dominance;
	solution_t  *solutions;
	ea_nsga2_t *nsga2_solutions_p;
	owner_file_t *solutions_ana = NULL; 
	char *column_1=NULL;
	char *column_2=NULL;
	int num_solutions;
	int num_obj;


	char *path_file_name;
	path_file_name = Malloc(char, MAX_PATH_FILE_NAME);
	strcpy(path_file_name, argv[1]);	
	column_1= Malloc(char, MAX_RANDOM_STRING);
	column_2= Malloc(char, MAX_RANDOM_STRING);

	//Getting Solutions
	solutions = loading_file_solutions(&num_solutions, &num_obj, path_file_name);
	in_param.size_population = num_solutions;
	in_param.number_fitness = num_obj;
	set_columns_of_solution_file(column_1, column_2, path_file_name);

/**************** START GETTING FRONT *************************/
	//Allocating solutions to analysis
	solutions_ana = loading_owner_file_solution(&num_solutions, &in_param.number_fitness, path_file_name);

	in_param.size_population = num_solutions;	
	nsga2_solutions_p = allocate_nsga2_without_allocation_of_representation(&in_param);	
	//Setting identification
	for (int i = 0; i < num_solutions; i++){
		nsga2_solutions_p[i].sol->ID = i+1;
	}

    //Setting dominance
	dominance = allocate_dominance(&num_solutions);
	set_dominance(dominance, solutions, &num_solutions); //solutions_p

	//Coping values of dominance
	for (int ind = 0; ind < num_solutions; ind++){
		// Indicates number of solutions that are dominated by me
		solutions_ana[ind].number_solutions_are_dominated = dominance[ind].max_dominated;
	}	


	//Coping values of objective
	for (int ind = 0; ind < num_solutions; ind++){
		for (int ob = 0; ob < in_param.number_fitness; ob++)	
			nsga2_solutions_p[ind].sol->obj_values[ob] = solutions[ind].obj_values[ob];
	}

    //Setting front based on dominance concept
    compute_fronts(nsga2_solutions_p, dominance, &num_solutions);

    //Coping values from nsga2_solutions_p to owner_file_t    
	for (int ind = 0; ind < num_solutions; ind++){
		/*
		//Coping objectives
		for (int ob = 0; ob < in_param.number_fitness; ob++){	
			solutions_ana[ind].obj_values[ob] = nsga2_solutions_p[ind].sol->obj_values[ob];
		}
		*/
		solutions_ana[ind].front = nsga2_solutions_p[ind].front;
	}
	desallocate_dominance(dominance, &num_solutions);
	desallocate_solution_nsga2(nsga2_solutions_p, &num_solutions);
/**************** FINISHED GETTING FRONT *************************/

/**************** START GETTING FINAL RESULTS *************************/
    //Sorting solutions
    sorting_solutions_by_front_dominance(solutions_ana, &num_solutions, &in_param.number_fitness);

    //Saving file
	save_analysis_files_no_objectives(solutions_ana, &num_solutions, &in_param.number_fitness, column_1, column_2);
/**************** FINISHED FINAL RESULTS *************************/
	free(column_1);
	free(column_2);
	desalocate_file_t(solutions_ana, &num_solutions);	
	desallocate_solution(solutions, &num_solutions);	
	free(path_file_name);

	display_msg("Done Computing Front !!! \n");

	return 0;
}