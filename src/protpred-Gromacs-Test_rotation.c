#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "defines.h"
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

int main(int argc, char *argv[]){
	input_parameters_t *in_para;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);

	FILE *f_angle = NULL;
    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    protein_t *population_p = NULL;
    float angle, angle_d;
    int num_res_first, chi, max_chi;
    int num_residue_choose_2;
    char *file_name, *file_name_aux;
    int num_frame, num_rotacoes, p, r, one;
    int kind_of_rotation;
    char *name_kid_of_rotation;
    protein_t *frame_rot; //frame for rotation
    protein_t *rotated_frames; //frames that were rotated

    //Allocating variables
    file_name = Malloc(char, 600);
    file_name_aux = Malloc(char, 600);
    name_kid_of_rotation = Malloc(char, 15);
    //Assing varible using for allocating 
    one = 1;
    num_frame = 1;

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating PDB ATOMS of population_p
    population_p = allocateProtein(&one);    

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &one, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);

    //Allocating frame
    frame_rot = allocateProtein(&num_frame);
    frame_rot->p_atoms = allocate_pdbatom(&population_p[0].p_topol->numatom);
    copy_protein(frame_rot, &population_p[0]);
    //Getting data to rotate
    printf ("Enter number of rotation: ");
    scanf ("%d",&num_rotacoes);        
    num_res_first = 0;
    while ( (num_res_first <= 0) || (num_res_first > population_p[0].p_topol->numres)){
        printf ("Enter number of residue for rotating: 1 to %d ", population_p[0].p_topol->numres);
        scanf ("%d",&num_res_first);            
    }
    angle_d = 0;
    while ( (angle_d <= 0) || (angle_d > 360)){
        printf ("Enter value of angle for rotating: 0 to 360 (degree) ");
        scanf ("%f",&angle_d);                    
    }    
    angle = degree2radians(&angle_d);
    kind_of_rotation = 0;
    while ( (kind_of_rotation <= 0) || (kind_of_rotation > 4)){
        printf ("Enter kid of rotating: 1 - PHI 2 - PSI 3 - OMEGA 4 - Side Chain ");
        scanf ("%d",&kind_of_rotation);                    
    }
    //Assing name of kind of rotation    
    if (kind_of_rotation == 1){
        strcpy(name_kid_of_rotation, "PHI");
    }else if (kind_of_rotation == 2){
        strcpy(name_kid_of_rotation, "PSI");
    }else if (kind_of_rotation == 3){
        strcpy(name_kid_of_rotation, "OMEGA");
    }if (kind_of_rotation == 4){
        strcpy(name_kid_of_rotation, "Side_Chain");
    }
    //Checking side chain rotation 
    if (kind_of_rotation == 4){
        char *res_name;
        res_name = Malloc(char, MAX_CHI);
        get_res_name_from_res_num(res_name, &num_res_first, population_p[0].p_atoms, &population_p[0].p_topol->numatom);        
        max_chi = get_number_chi(res_name);
        chi = 0;
        if (max_chi > 0){
            //Choose a chi of residue. It must be started 1 untill number of residue
             if (max_chi == 1){
                chi = 1;
             }else{                
                while ( (chi <= 0) || (chi > max_chi) ){                        
                    printf ("Enter Side Chain 1 to %d ", max_chi);
                    scanf ("%d",&chi);                    
                }
            }                    
        }else{
            printf("No side chain for your selected residue (Name %s Number %i). Please choose another residue\n",res_name, num_res_first);
            exit(1);
        }        
        free(res_name);
    }
    //Allocating frames for creating final file
    rotated_frames = allocateProtein(&num_rotacoes);
    for (r = 0; r < num_rotacoes; r++){
        rotated_frames[r].p_atoms = allocate_pdbatom(&population_p[0].p_topol->numatom);
    }                     
    //Performing the rotations
    for (p = 0; p < num_rotacoes; p++){
        if (kind_of_rotation == 1){
            rotation_phi(&frame_rot[0], &num_res_first, &angle);
        }else if (kind_of_rotation == 2){
            rotation_psi(&frame_rot[0], &num_res_first, &angle);
        }else if (kind_of_rotation == 3){
            rotation_omega(&frame_rot[0], &num_res_first, &angle);
        }else if (kind_of_rotation == 4){                        
            rotation_chi(&frame_rot[0], &num_res_first, &chi, &angle);
        }
        copy_protein(&rotated_frames[p], &frame_rot[0]);
    }
    //Saving rotated file 
    sprintf(file_name_aux,"2PG_Rotating_%i_Applying_%s_at_Res_%i.pdb",num_rotacoes,name_kid_of_rotation,num_res_first);
    strcpy(file_name, file_name_aux);
    save_population_file(rotated_frames, in_para->path_local_execute,
            file_name, &num_rotacoes );

    free(name_kid_of_rotation);
    desallocateProtein(population_p, &one);
    deAllocateload_parameters(in_para);

	return 0;
}
