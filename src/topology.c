#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defines.h"
#include "topology.h"
#include "consts.h"
#include "messages.h"
#include "topology_types.h"
#include "topologyio.h"
#include "topologylib.h"
#include "pdbatom.h"

#define TAM_BLOCO_PRM 80
#define NUM_HYDROGEN_BACKBONE_NITROGEN 4
#define NUM_PROTEIN_BACKBONE_ATOMS 8

/*Represents the Hydrogens of backbone that can be connected with N */
static const char *hydrogen_backbone_Nitrogen[] = {"H1", "H2", "H3", "HN"};

/*Represents the all atoms that can form the protein backbone*/
static const char *protein_backbone[] = {"N", "CA", "C", "O", "H1", "H2", "H3", "HN"};

static void desAllocate_top_residue_atom_info(top_residue_atom_info_t *info){	
	if (info->fixed_atoms != NULL){
		free(info->fixed_atoms);
	}
	if (info->moved_atoms != NULL){
		free(info->moved_atoms);	
	}
}

top_global_t *allocateTop_Global(const primary_seq_t *primary_sequence,
	const int *numatom){
	top_global_t *top_aux;

	top_aux = Malloc(top_global_t,1);
	top_aux->numatom = *numatom;
	top_aux->numres  = primary_sequence->num_res;

	top_aux->range_atoms  = Malloc(top_residue_range_atoms_t, top_aux->numres);
	top_aux->phi          = Malloc(top_residue_atom_info_t, top_aux->numres);
	top_aux->psi          = Malloc(top_residue_atom_info_t, top_aux->numres);
	top_aux->omega        = Malloc(top_residue_atom_info_t, top_aux->numres);
	top_aux->side_chains  = Malloc(top_residue_atom_info_t, top_aux->numres);

	return top_aux;
}
void  desAllocateTop_Global(top_global_t *top_aux){
	desAllocate_top_residue_atom_info(top_aux->phi);
	desAllocate_top_residue_atom_info(top_aux->psi);	
	desAllocate_top_residue_atom_info(top_aux->omega);
	free(top_aux->range_atoms);
	//desAllocate_top_residue_atom_info(top_aux->side_chains);
	free(top_aux);
}


static int get_atom_index_by_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname, const int *num_atom){
	const pdb_atom_t * aux;

	aux = search_pdb_atom_from_resnum_atom_name(atoms, res_num, atomname, num_atom);
	if (aux == NULL){
		fatal_error("Atom not found at get_atom_index_by_resnum_atom_name\n");
	}
	return aux->atmnumber;

}

static void build_topology_individual_omega(protein_t *prot){
	char *atom_N, *atom_C;
	int i_af;
	int next_res; // next residue

	atom_N = Malloc(char, 2);
	atom_C = Malloc(char, 2);	
	

	strcpy(atom_N, "N");
	strcpy(atom_C, "C");
	
	//The last residue (C-Terminal) does not make rotation
	for (int r = 1; r <= prot->p_topol->numres; r++){
		next_res = r + 1; //obtaing the next residue		
		if (next_res <= prot->p_topol->numres){
			//Build Fixed Atoms
			i_af = -1;
			//num_hydrogen_backbone = get_number_hydrogen_backbone(prot, &r);
			prot->p_topol->omega[r-1].num_fixed = 2;
			prot->p_topol->omega[r-1].fixed_atoms = Malloc(int, prot->p_topol->omega[r-1].num_fixed);
			i_af++;
			prot->p_topol->omega[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					&r, atom_C, &prot->p_topol->numatom);		
			i_af++;
			prot->p_topol->omega[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					&next_res, atom_N, &prot->p_topol->numatom);

			/*
			//Build Moved Atom		
			i_af = -1;
			//It is considered moved atom all atoms of NEXT residue (execept N because it is fixed)
			prot->p_topol->omega[r-1].num_moved = 2 + ((prot->p_topol->range_atoms[next_res-1].last_atom - prot->p_topol->range_atoms[next_res-1].first_atom) - prot->p_topol->omega[r-1].num_fixed);
			prot->p_topol->omega[r-1].moved_atoms = Malloc(int, prot->p_topol->omega[r-1].num_moved);
			for (int i_a = prot->p_topol->range_atoms[next_res-1].first_atom; i_a <= prot->p_topol->range_atoms[next_res-1].last_atom; i_a++){			
				if (is_fixed_atom(&prot->p_atoms[i_a-1].atmnumber, 
					prot->p_topol->omega[r-1].fixed_atoms, 
					&prot->p_topol->omega[r-1].num_fixed) == bfalse){				
					i_af++;				
					prot->p_topol->omega[r-1].moved_atoms[i_af] = prot->p_atoms[i_a-1].atmnumber;
				}
			
			}*/	
		}		
	}

	free(atom_N);	
	free(atom_C);	
}

static void build_topology_individual_psi(protein_t *prot){
	char *atom_C, *atom_CA, *atom_O;
	
	atom_CA = Malloc(char, 3);
	atom_C = Malloc(char, 2);	
	atom_O = Malloc(char, 2);
	

	strcpy(atom_CA, "CA");
	strcpy(atom_C, "C");
	strcpy(atom_O, "O");
	
	//The last residue does not make rotation
	for (int r = 1; r < prot->p_topol->numres; r++){
		//Build Fixed Atoms		
		prot->p_topol->psi[r-1].num_fixed = 2;
		prot->p_topol->psi[r-1].fixed_atoms = Malloc(int, prot->p_topol->psi[r-1].num_fixed);
		prot->p_topol->psi[r-1].fixed_atoms[0] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			&r, atom_CA, &prot->p_topol->numatom);
		prot->p_topol->psi[r-1].fixed_atoms[1] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			&r, atom_C, &prot->p_topol->numatom);
		//Build Moved Atom
		prot->p_topol->psi[r-1].num_moved = 1;
		prot->p_topol->psi[r-1].moved_atoms = Malloc(int, prot->p_topol->psi[r-1].num_moved);
		prot->p_topol->psi[r-1].moved_atoms[0] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			&r, atom_O, &prot->p_topol->numatom);
	}
	
	free(atom_CA);
	free(atom_C);
	free(atom_O);

}

static void build_topology_individual_phi(protein_t *prot){
	char *atom_N, *atom_C, *atom_CA, *atom_O;
	//int num_hydrogen_backbone;
	int i_af;

	atom_N = Malloc(char, 2);
	atom_CA = Malloc(char, 3);
	atom_C = Malloc(char, 2);	
	atom_O = Malloc(char, 2);

	strcpy(atom_N, "N");
	strcpy(atom_CA, "CA");
	strcpy(atom_C, "C");
	strcpy(atom_O, "O");

	//The first residue does not make rotation
	for (int r = 1; r <= prot->p_topol->numres; r++){
		if (r > 1){
			//Build Fixed Atoms
			i_af = -1;
			//num_hydrogen_backbone = get_number_hydrogen_backbone(prot, &r);
			prot->p_topol->phi[r-1].num_fixed = 2;//2 + num_hydrogen_backbone;
			prot->p_topol->phi[r-1].fixed_atoms = Malloc(int, prot->p_topol->phi[r-1].num_fixed);
			i_af++;
			prot->p_topol->phi[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					&r, atom_N, &prot->p_topol->numatom);
			i_af++;
			prot->p_topol->phi[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					&r, atom_CA, &prot->p_topol->numatom);
			//Build Moved Atom		
			i_af = -1;
			//It is considered moved atom all atoms that are not fixed atoms
			prot->p_topol->phi[r-1].num_moved = 1 + ((prot->p_topol->range_atoms[r-1].last_atom - prot->p_topol->range_atoms[r-1].first_atom) - prot->p_topol->phi[r-1].num_fixed);
			prot->p_topol->phi[r-1].moved_atoms = Malloc(int, prot->p_topol->phi[r-1].num_moved);
			for (int i_a = prot->p_topol->range_atoms[r-1].first_atom; i_a <= prot->p_topol->range_atoms[r-1].last_atom; i_a++){			
				if (is_fixed_atom(&prot->p_atoms[i_a-1].atmnumber, 
					prot->p_topol->phi[r-1].fixed_atoms, 
					&prot->p_topol->phi[r-1].num_fixed) == bfalse){				
					i_af++;				
					prot->p_topol->phi[r-1].moved_atoms[i_af] = prot->p_atoms[i_a-1].atmnumber;
				}
			}
		}else{
			prot->p_topol->phi[r-1].num_fixed = 0;
			prot->p_topol->phi[r-1].fixed_atoms = NULL;
			prot->p_topol->phi[r-1].num_moved = 0;
			prot->p_topol->phi[r-1].moved_atoms = NULL;
		}
	}

	free(atom_N);
	free(atom_CA);
	free(atom_C);
	free(atom_O);
}

static void build_topology_individual_side_chains(protein_t *prot){
	int i_af;
	for (int r = 1; r <= prot->p_topol->numres; r++){
		//Build Fixed Atoms - All backbone atoms
		i_af = -1;
		prot->p_topol->side_chains[r-1].num_fixed = get_number_atoms_backbone(prot, &r);
		prot->p_topol->side_chains[r-1].fixed_atoms = Malloc(int, prot->p_topol->side_chains[r-1].num_fixed);
		for (int i_a = prot->p_topol->range_atoms[r-1].first_atom; i_a <= prot->p_topol->range_atoms[r-1].last_atom; i_a++){
			if ( is_backbone_atom(prot->p_atoms[i_a-1].atmname) == btrue){
				i_af++;
				prot->p_topol->side_chains[r-1].fixed_atoms[i_af] = prot->p_atoms[i_a-1].atmnumber;
			}			
		}
		//Build Moved Atom - All atoms that are not backbone
		i_af = -1;
		prot->p_topol->side_chains[r-1].num_moved = 1 + ((prot->p_topol->range_atoms[r-1].last_atom - prot->p_topol->range_atoms[r-1].first_atom) - prot->p_topol->side_chains[r-1].num_fixed);
		prot->p_topol->side_chains[r-1].moved_atoms = Malloc(int, prot->p_topol->side_chains[r-1].num_moved);
		for (int i_a = prot->p_topol->range_atoms[r-1].first_atom; i_a <= prot->p_topol->range_atoms[r-1].last_atom; i_a++){
			if ( is_backbone_atom(prot->p_atoms[i_a-1].atmname) == bfalse){
				i_af++;
				prot->p_topol->side_chains[r-1].moved_atoms[i_af] = prot->p_atoms[i_a-1].atmnumber;
			}			
		}
	}
}

/** Builds the range of atoms per residue
*/
static void build_topology_individual_atoms(protein_t *prot){
	int residue_num_aux; // auxiliar for number of residue
	int i_a; // index atom
	int r; // index residue

	i_a = 0; 
	r = 0; 

	while (i_a < prot->p_topol->numatom){	
		residue_num_aux = prot->p_atoms[i_a].resnum;
		prot->p_topol->range_atoms[r].first_atom = prot->p_atoms[i_a].atmnumber; 
		// while is same residue do
		while (residue_num_aux == prot->p_atoms[i_a].resnum){ 
			prot->p_topol->range_atoms[r].last_atom = prot->p_atoms[i_a].atmnumber; 
			i_a = i_a + 1; // increase atom index
		}
		// Next residue
		r = r + 1; 
	}	
}

void build_topology_individual(protein_t *prot){
	build_topology_individual_atoms(prot);
	build_topology_individual_psi(prot);
	build_topology_individual_phi(prot);
	build_topology_individual_omega(prot);

	//build_topology_individual_side_chains(prot);
}

void build_topology_population(protein_t *pop, const int *pop_size){
	for (int i = 0; i < *pop_size; i++){
		build_topology_individual(&pop[i]);
	}
}

/** Returns the number of Hydrogens at backbone
*/
int get_number_hydrogen_backbone(const protein_t *prot, const int *numres){
	int num = 0;

	for (int i_a = prot->p_topol->range_atoms[*numres-1].first_atom; i_a <= prot->p_topol->range_atoms[*numres-1].last_atom; i_a++){
		if ( is_hydrogen_backbone_Nitrogen(prot->p_atoms[i_a-1].atmname) ){
			num = num + 1;
		}

	}
	return num;
}

/** Checks atom is a hydrogen at backbone that can be connected with N
*/
boolean_t is_hydrogen_backbone_Nitrogen(const char *atomname){	
		for (int i = 0; i < NUM_HYDROGEN_BACKBONE_NITROGEN; i++){
			if ( strncmp(atomname, hydrogen_backbone_Nitrogen[i], 2) == 0){
				return btrue;
			}
		}
		return bfalse;
}

/** Checks atmnumber is a fixed atom
*/
boolean_t is_fixed_atom(const int *atmnumber, const int *fixed_atoms, const int *num_fixed){
	for (int i = 0; i < *num_fixed; i++){
		if (*atmnumber == fixed_atoms[i]){
			return btrue;
		}
	}
	return bfalse;
}


/** Returns the number of atoms which are at backbone
*/
int get_number_atoms_backbone(const protein_t *prot, const int *numres){
	int num = 0;
	for (int i_a = prot->p_topol->range_atoms[*numres-1].first_atom; i_a <= prot->p_topol->range_atoms[*numres-1].last_atom; i_a++){
		if ( is_backbone_atom(prot->p_atoms[i_a-1].atmname) == btrue){
			num = num + 1;
		}

	}
	return num;
}


/** Checks atomname is a backbone atom
*/
boolean_t is_backbone_atom(const char *atomname){
	for (int i = 0; i < NUM_PROTEIN_BACKBONE_ATOMS ; i++){
		if ( strncmp(atomname, protein_backbone[i], 2) == 0){
			return btrue;
		}
	}
	return bfalse;
}
