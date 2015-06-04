#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
#include "diehdral.h"
#include "load_parameters.h"
#include "messages.h"
#include "defines.h"
#include "solution.h"
#include "futil.h"
#include "string_owner.h"
#include "protein.h"
#include "topology.h"
#include "pdbio.h"
#include "pdbatom.h"
#include "aminoacids.h"
#include "aminoacids_io.h"
#include "populationio.h"
#include "topologyio.h"
#include "rotation.h"
#include "math_owner.h"
#include "vector_math.h"



float radians2degree(const float *radians){
	return *radians*180/PI;
}

int main(int argc, char *argv[]){
	input_parameters_t *in_para;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);

	FILE *f_angle = NULL;
    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    protein_t *population_p = NULL;
    float *chi = NULL;
    int res_num, n_chi;
    float phi;
    float psi;
    float omega;    
    int chi_number;	    
    char file_name_angles[MAX_FILE_NAME];
    char *path_file_angle = NULL;

    int one = 1;

    //Allocanting chi for compute chi angles
    chi = (float*)malloc(sizeof(float)*MAX_CHI);

    //Seting file name angles
    strcpy(file_name_angles,"Dihedral_angles.txt");
    path_file_angle = path_join_file(in_para->path_local_execute, file_name_angles);
    f_angle = open_file(path_file_angle, fWRITE);

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating PDB ATOMS of population_p
    population_p = allocateProtein(&one);    

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &one, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);
	
    for (res_num = 1; res_num <= population_p[0].p_topol->numres; res_num++){    	
    	fprintf(f_angle,"phi\tpsi\tomega of res %d\n", res_num);
    	phi = compute_phi_residue(population_p[0].p_atoms, &res_num, population_p[0].p_topol);
    	psi = compute_psi_residue(population_p[0].p_atoms, &res_num, population_p[0].p_topol);
	   	omega = compute_omega_residue(population_p[0].p_atoms, &res_num, population_p[0].p_topol);
    	fprintf(f_angle,"%f\t%f\t%f\n", radians2degree(&phi), radians2degree(&psi), radians2degree(&omega) );
    	compute_chi_residue(chi, &chi_number, population_p[0].p_atoms, &res_num, population_p[0].p_topol);
		fprintf(f_angle,"Chi do residuo %d\n", res_num);	
		for (n_chi = 0; n_chi < chi_number; n_chi++){
			fprintf(f_angle,"%f\n", radians2degree(&chi[n_chi]) );
		}    	
	}
	fclose(f_angle);
    free(chi);
    free(path_file_angle);
    desallocateProtein(population_p, &one);
    deAllocateload_parameters(in_para);

	printf("Created %s file was success \n",file_name_angles);
	display_msg("Dihedral angles were computed with success \n");

	return 0;
}
