#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dominancep.h"
#include "dominancep_type.h"
#include "defines.h"
#include "messages.h"
#include "futil.h"


dominancep_t * allocate_dominancep(const int *size){
	dominancep_t *dominancep;
	dominancep = Malloc(dominancep_t, *size);
	return dominancep;
}

void desallocate_dominancep(dominancep_t *dominancep, const int *size){
	for (int d = 0; d < *size; d++){
		free(dominancep[d].set_dominated);
		free(dominancep[d].set_dominates_me);
	}
	free(dominancep);
	dominancep = NULL;
}

void show_dominancep(const dominancep_t *dominancep, const int *size){
	int i = 0;
	int j = 0;

	printf("\n");
	for (i = 0; i < *size; i++){
		printf("SOLUÇÃO %i\n", i+1);
		printf("How many solutions I dominate: %d\n",
			dominancep[i].amount_dominated);
		printf("What solutions I dominate:\t");
		for (j = 0; j < *size; j++){
			printf("%d\t\t", dominancep[i].set_dominated[j]);
		}
		printf("\n");
		printf("How many solutions dominate me: %d\n",
			dominancep[i].amount_dominates_me);
		printf("What solutions dominate me:\t");
		for(j = 0; j < *size; j++){
			printf("%d \t\t", dominancep[i].set_dominates_me[j]);
		}
		printf("\n");
		printf("Values of objectives:\t");
		for (j = 0; j < dominancep[i].sol->num_obj; j++){
			printf("obj %d: %.2f\t\t", j + 1, dominancep[i].sol->obj_values[j]);
		}
		printf("\n\n");

	}

}

void save_dominance(const dominancep_t *dominancep, const int *size){
	FILE *d_file=NULL;
	char *file_name;
	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "dominancep_output.txt");

	d_file = open_file(file_name, fWRITE);

	for (int d = 0; d < *size; d++){
/*
		fprintf(d_file, "Index dominance %i\n", d+1);
		fprintf(d_file, "How many solutions dominate me %i\n",
			dominancep[d].dominance.how_many_solutions_dominate_it);
		fprintf(d_file, "How many solutions are dominated by me %i\n",
			dominancep[d].dominance.max_dominated);
		fprintf(d_file, "Values of objective\n");
		for (int ob = 0; ob < dominancep[d].dominance.sol->num_obj; ob++){
			fprintf(d_file, "objective %i value is %f\n", ob,
				dominancep[d].dominance.sol->obj_values[ob]);
		}
		fprintf(d_file,"-------------------- \n\n");
*/
	}
	fclose(d_file);
	free(file_name);
}


/** Applies the dominancep concept in solutions
* dominancep is a pointer that stores the application of "robustez" concept in solution
* solutions array of Solutions
* size number of solutions
*/
void set_dominancep(dominancep_t *dominancep, const solution_t *solution, const int *size){

    initialize_dominancep(dominancep, solution, size);

    int dominated = 0;
    int dominates_me = 0;
    int i = 0;
    int j = 0;
    int o = 0;
    int qtde_objetivos = dominancep[0].sol->num_obj;

    for(i = 0; i < *size; i++){
	dominancep[i].set_dominated[i] = -1;
        dominancep[i].set_dominates_me[i] = -1;
        for(j = i + 1; j < *size; j++){
	    dominated = 0;
            dominates_me = 0;
	    //printf("\nSol %i compared to sol %i\n", i, j);
	    for(o = 0; o < qtde_objetivos; o++){
		//printf("Obj %i de i: %f\t", o + 1, dominancep[i].sol->obj_values[o]);
                //printf("Obj %i de j: %f\n", o + 1, dominancep[j].sol->obj_values[o]);
                if(dominancep[i].sol->obj_values[o] < dominancep[j].sol->obj_values[o]){
                    dominated++;
		    //printf("%i domina %i no objetivo %i.\n", i, j, o + 1);
                } else {
                    if(dominancep[i].sol->obj_values[o] > dominancep[j].sol->obj_values[o]){
                        dominates_me++;
		        //printf("%i é dominado por %i no objetivo %i.\n", i, j, o + 1);

                    }
                }
            }
	    //printf("dominated = %d\t\tdominates_me = %d\n\n", dominated, dominates_me);
            if((dominated >= 1) && (dominates_me == 0)){
                //Talvez seja necessário alterar essa implementação a fim de
                //trabalhar com populacão externa do spea2.
                //Nesse caso penso em utilizar o ID de solution_t;
		//printf("%i domina %i.\n", i, j);
                dominancep[i].set_dominated[j] = 1;
		dominancep[j].set_dominated[i] = -1;
		dominancep[i].amount_dominated++;

                dominancep[j].set_dominates_me[i] = 1;
                dominancep[i].set_dominates_me[j] = -1;
		dominancep[j].amount_dominates_me++;

            } else {
                if(dominates_me >= 1 && dominated == 0){
			//printf("%i é dominado por %i.\n", i, j);
			dominancep[j].set_dominated[i] = 1;
			dominancep[i].set_dominated[j] = -1;
			dominancep[j].amount_dominated++;

	            	dominancep[i].set_dominates_me[j] = 1;
	            	dominancep[j].set_dominates_me[i] = -1;
			dominancep[i].amount_dominates_me++;
                } else {
		        //printf("ninguém domina ninguém entre %i e %i.\n", i, j);
			dominancep[i].set_dominated[j] = -1;
			dominancep[j].set_dominated[i] = -1;
			dominancep[i].set_dominates_me[j] = -1;
			dominancep[j].set_dominates_me[i] = -1;
		    }

                }
            }
	}
}

/** Initialize dominancep struct
* dominancep is a pointer that stores the application of dominancep concept in solution
* solutions array of Solutions
* size number of solutions
*/
static void initialize_dominancep(dominancep_t *dominancep, const solution_t *solutions, const int *size){
    int s = 0;
    for(s = 0; s < *size; s++){
		dominancep[s].sol = &solutions[s];
		dominancep[s].set_dominated = Malloc(int, *size);
		dominancep[s].set_dominates_me = Malloc(int, *size);
		dominancep[s].amount_dominated = 0;
		dominancep[s].amount_dominates_me = 0;
	}
}


