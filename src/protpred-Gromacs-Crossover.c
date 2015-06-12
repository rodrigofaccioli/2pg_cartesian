#include <stdlib.h>
#include <string.h>
#include <math.h>

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
#include "randomlib.h"
#include "diehdral.h"


void crossover_diehdral(protein_t *new_prot, const protein_t * p1, const protein_t * p2) {
    int res_num_choose, next_res, max_res, res, a, chi_ref, num_chi;   
    float phi_1, phi_2, psi_1, psi_2, omega_1, omega_2, chi_1, chi_2,angle_diff_phi, angle_diff_psi, angle_diff_omega, angle_diff_chi;    
    float *v_chi_1, *v_chi_2;
    int side_chain;
    max_res = get_number_res_from_atom(new_prot->p_atoms, &new_prot->p_topol->numatom);  
    res_num_choose = 10;//_get_int_random_number(&max_res);

    printf("%d\n", res_num_choose);

    //Allocanting chi for compute chi angles
    v_chi_1 = (float*)malloc(sizeof(float)*MAX_CHI);
    v_chi_2 = (float*)malloc(sizeof(float)*MAX_CHI);


    //Coping atoms from p1 to new_prot
    a = 0;
    for (res = p1->p_atoms[a].resnum; res <= max_res; res++){        
        //Coping all atoms of residue res
        while (p1->p_atoms[a].resnum == res){
            copy_pdb_atom(&new_prot->p_atoms[a], &p1->p_atoms[a]);
            a = a + 1;                
        }
    }
    //Next residue of choose residue
    next_res = res_num_choose + 1;
    //Performing the rotation in remaining protein 
    for (res = next_res; res <= max_res; res++){
        //Rotating residue
        //phi
        phi_1 = compute_phi_residue(p1->p_atoms, &res, p1->p_topol);
        phi_2 = compute_phi_residue(p2->p_atoms, &res, p2->p_topol);
        angle_diff_phi = phi_2 - phi_1;
        rotation_phi(new_prot, &res, &angle_diff_phi); 
        //psi
        psi_1 = compute_psi_residue(p1->p_atoms, &res, p1->p_topol);
        psi_2 = compute_psi_residue(p2->p_atoms, &res, p2->p_topol);
        angle_diff_psi = psi_2 - psi_1;
        rotation_psi(new_prot, &res, &angle_diff_psi);
        //omega
        omega_1 = compute_omega_residue(p1->p_atoms, &res, p1->p_topol);
        omega_2 = compute_omega_residue(p2->p_atoms, &res, p2->p_topol);
        angle_diff_omega = omega_2 - omega_1;
        rotation_omega(new_prot, &res, &angle_diff_omega);
        //chi        
        compute_chi_residue(v_chi_1, &num_chi, p1->p_atoms, &res, p1->p_topol);
        compute_chi_residue(v_chi_2, &num_chi, p2->p_atoms, &res, p2->p_topol);
        if (num_chi > 0){
            for (chi_ref = 0; chi_ref < num_chi; chi_ref++){
                side_chain = chi_ref + 1;
                angle_diff_chi = v_chi_2[chi_ref] - v_chi_1[chi_ref];
                rotation_chi(new_prot, &res, &side_chain, &angle_diff_chi);
            }
        }
    }
    free(v_chi_1);
    free(v_chi_2);
}

int main(int argc, char *argv[]){
	input_parameters_t *in_para;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);

    char *crossover_file_name = NULL;
    int ind, ind_ref_1, ind_ref_2;

    //Setting crossover file name 
    crossover_file_name = Malloc(char, MAX_FILE_NAME);
    strcpy(crossover_file_name, "pop_crossover.pdb");

    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    
    protein_t *population_p; // main population
    protein_t *new_population; // main population

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating PDB ATOMS of population_p
    population_p = allocateProtein(&in_para->size_population);
    new_population  = allocateProtein(&in_para->size_population);   

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);
    //Only to initialize new_population
    load_initial_population_file(new_population, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);
    initialize_protein_population_atoms(new_population, &in_para->size_population);

    for (ind = 0; ind < in_para->size_population; ind++){
        ind_ref_1 = 0;
        ind_ref_2 = 1;
        crossover_diehdral(&new_population[ind], &population_p[ind_ref_1], &population_p[ind_ref_2]);
    }

    save_population_file(new_population, in_para->path_local_execute, crossover_file_name, &in_para->size_population);

	deAllocateload_parameters(in_para);

	display_msg("Done 2PG Crossover \n");
	return 0;
}