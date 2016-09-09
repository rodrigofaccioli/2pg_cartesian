#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "owner_file_analysis.h"
#include "analysis_types.h"
#include "dominancep.h"
#include "dominancep_type.h"
#include "nsga2.h"
#include "defines.h"
#include "messages.h"
#include "futil.h"


#define ARRED(x) floor( x + 0.5 ) //Isso retorna um inteiro
#define ARRED2(x, z) floor( x*pow(10, z) + 0.5 )/pow(10, z) //sendo z o numero de casas decimais

static int compare_front_nsga2( const void *x, const void *y){
    nsga2_t * fx = (nsga2_t*) x;
    nsga2_t * fy = (nsga2_t*) y;

    if (fx->front > fy->front){
        return 1;
    } else { 
	if(fx->front == fy->front) {
	    if(fx->dominancep->amount_dominated < fy->dominancep->amount_dominated){
		return 1;
	    } else {
		return 0;
	    }
	} else {
	    return 0;	
	}    
    }
}

/*
static int compare_front_nsga2( const nsga2_t * nsga2, const nsga2_t * nsga2){
    if (nsga2->front > nsga2->front){
        return 1;
    } else { 
	if(nsga2->front == nsga2->front) {
	    if(nsga2->dominancep->amount_dominated < nsga2->dominancep->amount_dominated){
		return 1;
	    }
	} else {
	    return 0;	
	}    
    }
}
*/
nsga2_t * allocate_nsga2(const int *size){
	nsga2_t *nsga2;
	nsga2 = Malloc(nsga2_t, *size);
	return nsga2;
}

void desallocate_nsga2(nsga2_t *nsga2){
	free(nsga2);
	nsga2 = NULL;
}

void show_nsga2(const nsga2_t *nsga2, const int *size){
/*	
	int i = 0;
	int j = 0;
	int o = 0;

	printf("\n");
	for (i = 0; i < *size; i++){				
		printf("SOLUÃÃO %i\n", i+1);
		printf("How many solutions I dominate: %i\n",
			dominancep[i].amount_dominated);
		printf("What solutions I dominate:\t");
		for (j = 0; j < *size; j++){
			printf("%i\t", dominancep[i].set_dominated[j]);
		}
		printf("\n");
		printf("How many solutions dominate me %i\n",
			dominancep[i].amount_dominates_me);
		printf("What solutions dominate me:\n");
		for(j = 0; j < *size; j++){
			printf("Sol %i\t", dominancep[i].set_dominates_me[j]);
		}
		printf("Values of objectives\n");
-		for (o = 0; o < dominancep[i].sol->num_obj; o++){
			printf("obj %i: %f\t", o + 1,	dominancep[i].sol->obj_values[o]);
		}
		printf("\n\n");

	}
*/
}

void save_nsga2_v1(const nsga2_t *nsga2, const int *size){

	FILE *d_file=NULL;
	char *file_name;
	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "nsga2_output_in_column.txt");
	d_file = open_file(file_name, fWRITE);
	fprintf(d_file, "#PDB files of protein struct prediction ranked\n");
	fprintf(d_file, "#according to NSGA2 mechanism of ordering.\n");
	fprintf(d_file, "#This means that a multiobjective optimization\n");		
	fprintf(d_file, "#criteria was used\n\n");
	fprintf(d_file, "Rank\t");
	fprintf(d_file, "ID\t");
	fprintf(d_file, "File Name \t\t");
	fprintf(d_file, "Dominated \t");
	fprintf(d_file, "Dominates \t");
	fprintf(d_file, "Front\t\t");
	fprintf(d_file, "\n");
	fprintf(d_file, "\t\t\t\t\tby me\t");
	fprintf(d_file, "\tme\t");

	fprintf(d_file, "\n\n");
	for (int i = 0; i < *size; i++){
		fprintf(d_file, "%d\t", i+1);
		fprintf(d_file, "%d\t", nsga2[i].dominancep->sol->ID);
		fprintf(d_file, "%.*s\t\t", 15, nsga2[i].file_name);
		fprintf(d_file, "%d\t\t", nsga2[i].dominancep->amount_dominated);
		fprintf(d_file, "%d\t\t", nsga2[i].dominancep->amount_dominates_me);
		fprintf(d_file, "%d\t\t", nsga2[i].front);
		fprintf(d_file, "\n");
	}

	fclose(d_file);
	free(file_name);
}

void save_nsga2_v2(const nsga2_t *nsga2, const int *size){
/*
	FILE *d_file=NULL;
	char *file_name;
	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "nsga2_output.txt");

	d_file = open_file(file_name, fWRITE);
	
	for (int i = 0; i < *size; i++){

		fprintf(d_file, "Solution ID: %d\n", nsga2[i].dominancep->sol->ID);

		fprintf(d_file, "How many solutions are dominated by me %i\n",
			nsga2[i].dominancep->amount_dominated);
		fprintf(d_file, "How many solutions dominate me %i\n",
			nsga2[i].dominancep->amount_dominates_me);
		fprintf(d_file, "Fitness value: %f\n", nsga2[i].fitness);
		fprintf(d_file, "Value of objectives\n");
		for (int ob = 0; ob < nsga2[i].dominancep->sol->num_obj; ob++){
			fprintf(d_file, "objective %i value is %f\n", ob + 1,
				nsga2[i].dominancep->sol->obj_values[ob]);
		}
		fprintf(d_file,"\n-------------------- \n\n");

	}
	fclose(d_file);
	free(file_name);
*/
}

/** Calculates raw_fitness for each solution
* nsga2 is a pointer that stores data to calculate de "robustez" of a solution
* dominancep is a pointer that carries the values of dominance of a solution
* size number of solutions
*/
void set_front(nsga2_t *nsga2, const dominancep_t *dominancep, const int *size){

    
    int i = 0;
    int j = 0;
    int k = 0;
    int begin = 0;
    int end = 0;
    int front_aux = 0;
    int pos = 0;
    int pos_aux = 0;
    int done_size = 0;
    int aux = 0;
    int count_zero = 0;
    int count_one = 0;
    int done[*size];
    int mat[*size][*size];
    
    for(i = 0; i < *size; i++){
	for( j = 0; j < *size; j++){
	    if(nsga2[i].dominancep->set_dominates_me[j] == 1){
		mat[i][j] = 1;
	    } else {
	        mat[i][j] = 0;
            }
	}
    }

   while(k < *size){
	aux = 0;
        for(i = 0; i < *size; i++){
	    count_zero = 0;
	    count_one = 0;
	    for(j = 0; j < *size; j++){
	        if(mat[i][j] == 1){
		    count_one++;
		} else {
		    if(mat[i][j] == 0){
		        count_zero++;	  
		    }
		}
	    }
	    if(count_one == 0 && count_zero > 0){
		nsga2[i].front = front_aux;
		aux++;
		done[pos_aux++] = i;
	    }	
	}
	front_aux++;
	begin = end;
	end += aux;
	for(i = begin; i < end; i++){
	    for(j = 0; j < *size; j++){
		pos = done[i];
	        mat[pos][j] = -1;
	        mat[j][pos] = -1;
	    }	        
	}

	k += aux;
    }
}

/** Calculates the density for each solution
* nsga2 is a pointer that stores data to calculate the "robustez" of a solution
* dominancep is a pointer that carries the values of dominance of a solution
* size is the number of solutions
* size_external is the number of solutions on the external population file
*/



void order(nsga2_t *nsga2, const int *size){
    qsort(nsga2, *size, sizeof(nsga2_t), &compare_front_nsga2);
    //for(int i = 0; i < *size; i++){
	//printf("ID: %d -> front: %d\n", nsga2[i].dominancep->sol->ID, nsga2[i].front);
	//printf("mat[%d] = %f", i, );
    //}
    
    printf("\n");
    
}

/** Initialize dominancep struct
* dominancep is a pointer that stores the application of dominancep concept in solution
* solutions array of Solutions
* size number of solutions
*/
void initialize_nsga2(nsga2_t *nsga2, const dominancep_t *dominancep, const int *size, const owner_file_t * file_names){
    int s = 0;
    for(s = 0; s < *size; s++){
		nsga2[s].file_name = file_names[s].file_name;
		nsga2[s].dominancep = &dominancep[s];
		nsga2[s].front = -1;
		nsga2[s].crowding_distance = 0;
		nsga2[s].ranking = 0;
	}
}


