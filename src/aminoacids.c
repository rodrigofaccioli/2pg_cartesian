#include <string.h>

#include "aminoacids.h"

static void initialize_primary_seq(primary_seq_t* seq){

	for (int i = 0; i < seq->num_res; i++){
		strcpy(seq->seq_res[i].id_1, "");
		strcpy(seq->seq_res[i].id_3, "");
		seq->id = aNR;		
	}

}

primary_seq_t* allocate_primary_seq(const int *num_res){
	primary_seq_t* aux;

	aux = Malloc(primary_seq_t, 1);
	aux->num_res = *num_res;
	aux->seq_res = Malloc(amino_t, aux->num_res);

	initialize_primary_seq(aux);

	return aux;
}

void desallocate_primary_seq(primary_seq_t* seq){
	if (seq->seq_res != NULL){
		free(seq->seq_res);
	}
	free(seq);
}