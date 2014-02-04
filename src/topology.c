#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defines.h"
#include "enums.h"
#include "protein.h"
#include "topology.h"
#include "consts.h"
#include "messages.h"
#include "topology_types.h"
#include "topologyio.h"
#include "topologylib.h"

#define TAM_BLOCO_PRM 80


top_global_t *allocateTop_Global(const primary_seq_t *primary_sequence,
	const int *numatom){
	top_global_t *top_aux;

	top_aux = Malloc(top_global_t,1);
	top_aux->numatom = *numatom;
	top_aux->numres = primary_sequence->num_res;
	//top_aux->top_global_atom = Malloc(top_global_atom_t,*numatom);
	top_aux->top_global_dieh_phi = Malloc(top_global_dihedral_t,top_aux->numres);
	top_aux->top_global_dieh_psi = Malloc(top_global_dihedral_t,top_aux->numres);
	top_aux->top_global_dieh_omega = Malloc(top_global_dihedral_t,top_aux->numres);	
	//top_aux->top_global_res_atm = Malloc(top_global_res_atm_t,*numres);
	//top_aux->top_global_res_atms_bond = Malloc(top_global_res_atms_bond_t,*numatom);
	//top_aux->top_global_res_atms_bond_angle = Malloc(top_global_res_atms_bond_angle_t,*num_bond_angles);
	//top_aux->top_global_dieh_side_chains = Malloc(top_global_dihedral_side_chain_t,*num_side_chains);
	//top_aux->top_global_dihedral_angles_type = Malloc(top_global_dihedral_angles_type_t, *number_dihedrals_type);	
	return top_aux;
}
void  desAllocateTop_Global(top_global_t *top_global){
	//Falta criar um correto desallocate	
	free(top_global->top_global_dieh_phi);
	free(top_global->top_global_dieh_psi);
	free(top_global->top_global_dieh_omega);
/*
	free(top_global->top_global_atom);
	free(top_global->top_global_res_atm);
	free(top_global->top_global_res_atms_bond);
	free(top_global->top_global_res_atms_bond_angle);
	free(top_global->top_global_dieh_side_chains);
	free(top_global->top_global_dihedral_angles_type);
*/	
	free(top_global);
}


//void build_topology_pdb_population(pdb_atom_t** pop, const int *pop_size, 
//	const primary_seq_t *primary_sequence, const int *num_atom){

//}
