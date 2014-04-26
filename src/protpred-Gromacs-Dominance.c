#include <stdlib.h>
#include <string.h>

#include "dominance.h"
#include "load_parameters.h"
#include "messages.h"
#include "defines.h"
#include "solution.h"
#include "futil.h"
#include "string_owner.h"

static dominance_t *dominance;
static solution_t  *solutions;
static int num_solutions;
static int num_obj;

static void initialize_for_reading_file(const char *path_file_name){
	FILE * dominance_file;	
	char *line;
	char *line_splited;
	

	num_solutions = 0;
	num_obj = -1; //Because the first collum is index.
	line = Malloc(char, MAX_LINE_DOMINANCE);
	dominance_file = open_file(path_file_name, fREAD);

	//getting first line to obtain number of objectivies
	fgets(line,MAX_LINE_DOMINANCE,dominance_file);
 	line_splited = strtok (line,"\t");
  	while (line_splited != NULL){
    	num_obj = num_obj + 1;
    	line_splited = strtok(NULL, "\t");
  	}	
	//Counting how many lines the file has
	while ( fgets(line,MAX_LINE_DOMINANCE,dominance_file) != NULL){
		num_solutions = num_solutions + 1;
	}

	fclose(dominance_file);
	free(line);
}

static void loading_file_dominance(const char *path_file_name){
	FILE * dominance_file;	
	char *line;
	char *line_splited;
	int sol;

	line = Malloc(char, MAX_LINE_DOMINANCE);

	//Getting initial information
	initialize_for_reading_file(path_file_name);

	//Alocating Solution
	solutions = allocate_solution(&num_solutions, &num_obj);

	//Reading file and set values of objective
	sol = -1;
	dominance_file = open_file(path_file_name, fREAD);
	//Removing first line that is collumn
	fgets(line,MAX_LINE_DOMINANCE,dominance_file);
	while ( fgets(line,MAX_LINE_DOMINANCE,dominance_file) != NULL){
		sol = sol + 1;
		//Removing index collumn
		line_splited = strtok (line,"\t");
		//Setting number of objectives
		solutions[sol].num_obj = num_obj;
		for (int ob = 0; ob < num_obj; ob++){
			line_splited = strtok (NULL,"\t");
			trim(line_splited);
			solutions[sol].obj_values[ob] = str2double(line_splited);
		}					
	}	
	fclose(dominance_file);

	free(line);
}

int main(int argc, char *argv[]){
	input_parameters_t in_param;

	char *path_file_name;
	path_file_name = Malloc(char, MAX_PATH_FILE_NAME);
	strcpy(path_file_name, argv[1]);	

	loading_file_dominance(path_file_name);

	dominance = allocate_dominance(&num_solutions);
	set_dominance(dominance, solutions, &num_solutions);
	show_dominance(dominance, &num_solutions);
	
	desallocate_dominance(dominance, &num_solutions);
	desallocate_solution(solutions, &num_solutions);
	free(path_file_name);

	display_msg("Done Dominance !!! \n");

	return 0;
}