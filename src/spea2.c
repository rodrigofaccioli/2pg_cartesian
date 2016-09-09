#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//#include "owner_file_analysis.h"
#include "analysis_types.h"
#include "dominancep.h"
#include "dominance.h"
#include "dominance_type.h"
#include "spea2.h"
#include "defines.h"
#include "messages.h"
#include "futil.h"


#define ARRED(x) floor( x + 0.5 ) //Isso retorna um inteiro
#define ARRED2(x, z) floor( x*pow(10, z) + 0.5 )/pow(10, z) //sendo z o numero de casas decimais

static int compare_fitness( const void *x, const void *y){
    spea2_t * fx = (spea2_t*) x;
    spea2_t * fy = (spea2_t*) y;

    if (fx->fitness > fy->fitness){
        return 1;
    }else {
        return 0;
    }
}

int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

spea2_t * allocate_spea2(const int *size){
	spea2_t *spea2;
	spea2 = Malloc(spea2_t, *size);
	return spea2;
}

void desallocate_spea2(spea2_t *spea2){
	free(spea2);
	spea2 = NULL;
}

void show_spea2(const spea2_t *spea2, const int *size){
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
		for (o = 0; o < dominancep[i].sol->num_obj; o++){
			printf("obj %i: %f\t", o + 1,	dominancep[i].sol->obj_values[o]);
		}
		printf("\n\n");

	}
*/
}

void save_spea2_v1(const spea2_t *spea2, const int *size){

	FILE *d_file=NULL;
	char *file_name;
	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "spea2_output_in_column.txt");
	d_file = open_file(file_name, fWRITE);
	fprintf(d_file, "#PDB files of protein struct prediction ranked\n");
	fprintf(d_file, "#according toSPEA2 mechanism of ordering.\n");
	fprintf(d_file, "#This means that a multiobjective optimization\n");		
	fprintf(d_file, "#criteria was used\n\n");
	fprintf(d_file, "Rank\t");
	fprintf(d_file, "ID\t");
	fprintf(d_file, "File Name \t\t");
	fprintf(d_file, "Dominated \t");
	fprintf(d_file, "Dominates \t");
	fprintf(d_file, "Raw fit\t\t");
	fprintf(d_file, "Density\t\t");
	fprintf(d_file, "Fitness");
	fprintf(d_file, "\n");
	fprintf(d_file, "\t\t\t\t\tby me\t");
	fprintf(d_file, "\tme\t");
	/*
	for (int ob = 0; ob < spea2[0].dominancep->sol->num_obj; ob++){
		fprintf(d_file, "objective %d\t", ob + 1);
	}
	*/
	fprintf(d_file, "\n\n");
	for (int i = 0; i < *size; i++){
		fprintf(d_file, "%d\t", spea2[i].ranking);
		fprintf(d_file, "%d\t", spea2[i].dominancep->sol->ID);
		fprintf(d_file, "%.*s\t\t", 15, spea2[i].file_name);
		fprintf(d_file, "%d\t\t", spea2[i].dominancep->amount_dominated);
		fprintf(d_file, "%d\t\t", spea2[i].dominancep->amount_dominates_me);
		fprintf(d_file, "%.5f\t\t", spea2[i].raw_fitness);
		fprintf(d_file, "%.5f\t\t", spea2[i].density);
		fprintf(d_file, "%.5f", spea2[i].fitness);
/*		
		for (int ob = 0; ob < spea2[i].dominancep->sol->num_obj; ob++){
			fprintf(d_file, "%f\t",	spea2[i].dominancep->sol->obj_values[ob]);
		}	
*/
		fprintf(d_file, "\n");
	}

	fclose(d_file);
	free(file_name);

}

void save_spea2_v2(const spea2_t *spea2, const int *size){

	FILE *d_file=NULL;
	char *file_name;
	file_name = Malloc(char, MAX_FILE_NAME);
	strcpy(file_name, "spea2_output.txt");

	d_file = open_file(file_name, fWRITE);
	
	for (int i = 0; i < *size; i++){

		fprintf(d_file, "Solution ID: %d\tFile_name: %s\n", spea2[i].dominancep->sol->ID, spea2[i].file_name);
		fprintf(d_file, "How many solutions are dominated by me %i\n", spea2[i].dominancep->amount_dominated);
		fprintf(d_file, "How many solutions dominate me %i\n", spea2[i].dominancep->amount_dominates_me);
		fprintf(d_file, "Fitness value: %f\n", spea2[i].fitness);
		fprintf(d_file, "Value of objectives\n");
		for (int ob = 0; ob < spea2[i].dominancep->sol->num_obj; ob++){
			fprintf(d_file, "objective %i value is %f\n", ob + 1,
				spea2[i].dominancep->sol->obj_values[ob]);
		}
		fprintf(d_file,"\n-------------------- \n\n");

	}
	fclose(d_file);
	free(file_name);

}

/** Calculates raw_fitness for each solution
* spea2 is a pointer that stores data to calculate de "robustez" of a solution
* dominancep is a pointer that carries the values of dominance of a solution
* size number of solutions
*/
void set_raw_fitness(spea2_t *spea2, const dominancep_t *dominancep, const int *size){

    
    int i = 0;
    int j = 0;
    int k = 0;

    
    for(i = 0; i < *size; i++){
	for( j = 0; j < *size; j++){
	    if(spea2[i].dominancep->set_dominates_me[j] == 1){
		spea2[i].raw_fitness += spea2[j].dominancep->amount_dominated;
	    }
	}
	//printf("spea2[%d].raw_fitness = %f\n", i, spea2[i].raw_fitness);
    }
    //printf("\n");
}

/** Calculates the density for each solution
* spea2 is a pointer that stores data to calculate the "robustez" of a solution
* dominancep is a pointer that carries the values of dominance of a solution
* size is the number of solutions
* size_external is the number of solutions on the external population file
*/

void set_density(spea2_t *spea2, const int *size, const int *size_external){
    int i = 0;
    int j = 0;
    int o = 0;
    int obj = spea2[0].dominancep->sol->num_obj;
    double aux_soma = 0;
    double distance = 0;
    double k_element = 0;
    int k_element_int = 0;
    double mat[*size][*size];

    k_element = ARRED(sqrt(*size + *size_external));
    k_element_int = (int)k_element - 1;
    //printf("k_element = %.0f.\n", k_element);
	
    for(i = 0; i < *size; i++){
	mat[i][i] = 0;
	for( j = i+1; j < *size; j++){
	    aux_soma = 0;
	    for(o = 0; o < obj; o++){
		aux_soma += pow((spea2[i].dominancep->sol->obj_values[o] - spea2[j].dominancep->sol->obj_values[o]),2);	
	    }
	    //printf("aux_soma = %f\n", aux_soma);
	    distance = sqrt(aux_soma);
	    //printf("distance entre %d e %d = %f\n", i, j, distance); 
	    mat[i][j] = distance;
	    //printf("matriz[%d][%d] = %f\t\t", i, j, mat[i][j]);    
	    mat[j][i] = distance;	
	    //printf("matriz[%d][%d] = %f\t\n\n", j, i, mat[i][j]);    
	}
	qsort(mat[i], *size, sizeof(double), compare);
	//printf("####matriz[%d][%d] = %.2f\t\n", i, k_element_int, mat[i][(int)k_element]);
	spea2[i].density = 1/(mat[i][k_element_int] + 2);
	//printf("spea2[%d].density = %.30f\n", i, spea2[i].density);
    }
    printf("\n");
}


void set_fitness(spea2_t *spea2, const int *size){

    for(int i = 0; i < *size; i++){
	spea2[i].fitness = spea2[i].density + spea2[i].raw_fitness;
	//printf("spea2[%d].fitness = %0.30f\n", i, spea2[i].fitness); 
    }
    printf("\n");
}

void order(spea2_t *spea2, const int *size){

    /*
    double mat[*size];
    for(int i = 0; i < *size; i++){
	mat[i] = spea2[i].fitness;
    }
    qsort(mat, *size, sizeof(double), compare);
    */
    qsort(spea2, *size, sizeof(spea2_t), compare_fitness);
    for(int i = 0; i < *size; i++){
	spea2[i].ranking = i+1;
    }
    
}

/** Initialize dominancep struct
* dominancep is a pointer that stores the application of dominancep concept in solution
* solutions array of Solutions
* size number of solutions
*/
void initialize_spea2(spea2_t *spea2, const dominancep_t *dominancep, const int *size, const owner_file_t * file_names){
    int s = 0;
    for(s = 0; s < *size; s++){
		spea2[s].file_name = file_names[s].file_name;
		spea2[s].dominancep = &dominancep[s];
		spea2[s].raw_fitness = 0;
		spea2[s].density = 0;
		spea2[s].fitness = 0;
		spea2[s].ranking = 0;
	}
}


