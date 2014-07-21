#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "defines.h"
#include "enums.h"
#include "topology.h"
#include "topologylib.h"
#include "consts.h"
#include "messages.h"
#include "topology_types.h"
#include "topologyio.h"
#include "topologylib.h"
#include "pdbatom.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

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
		info->fixed_atoms = NULL;
	}
	if (info->moved_atoms != NULL){
		free(info->moved_atoms);	
		info->moved_atoms = NULL;
	}
}

static void desAllocate_top_residue_side_chains(top_residue_side_chains_t *side_chains,
	const int *numres){
	for (int r = 0; r <*numres;r++){
		if (side_chains[r].num_chi > 0){
			desAllocate_top_residue_atom_info(side_chains[r].atoms_chi);
		}
	}

}

_2PG_CARTESIAN_EXPORT
top_global_t *allocateTop_Global(const int *numres,
	const int *numatom){
	top_global_t *top_aux;

	top_aux = Malloc(top_global_t,1);
	top_aux->numatom = *numatom;
	top_aux->numres  = *numres; //primary_sequence->num_res;

	top_aux->range_atoms  = Malloc(top_residue_range_atoms_t, top_aux->numres);
	top_aux->phi          = Malloc(top_residue_atom_info_t, top_aux->numres);
	top_aux->psi          = Malloc(top_residue_atom_info_t, top_aux->numres);
	top_aux->omega        = Malloc(top_residue_atom_info_t, top_aux->numres);
	top_aux->side_chains  = Malloc(top_residue_side_chains_t, top_aux->numres);

	return top_aux;
}

_2PG_CARTESIAN_EXPORT
void  desAllocateTop_Global(top_global_t *top_aux){
	desAllocate_top_residue_atom_info(top_aux->phi);
	desAllocate_top_residue_atom_info(top_aux->psi);	
	desAllocate_top_residue_atom_info(top_aux->omega);
	desAllocate_top_residue_side_chains(top_aux->side_chains, &top_aux->numres);
	free(top_aux->range_atoms);	
	free(top_aux);
	top_aux = NULL;
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
			prot->p_topol->omega[r-1].num_fixed = 2;
			prot->p_topol->omega[r-1].fixed_atoms = Malloc(int, prot->p_topol->omega[r-1].num_fixed);
			i_af++;
			prot->p_topol->omega[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					&r, atom_C, &prot->p_topol->numatom);		
			i_af++;
			prot->p_topol->omega[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					&next_res, atom_N, &prot->p_topol->numatom);
			
			prot->p_topol->omega[r-1].num_moved = 0;
			prot->p_topol->omega[r-1].moved_atoms = NULL;

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
	int do_residue;

	atom_CA = Malloc(char, 3);
	atom_C = Malloc(char, 2);	
	atom_O = Malloc(char, 2);
	

	strcpy(atom_CA, "CA");
	strcpy(atom_C, "C");
	strcpy(atom_O, "O");
	
	//The last residue does not make rotation
	for (int r = 1; r < prot->p_topol->numres; r++){
		do_residue = 1;
		if (r == 1){
			/* When residue is N-Termninal have is necessary check 
			* residue is ACE or NME. These residues are special.  
			*/
			if ( (strcmp(prot->p_atoms[0].resname, "ACE") == 0) ||
				 (strcmp(prot->p_atoms[0].resname, "NME") == 0) ){
				//Build Fixed and moved Atoms for ACE and NME
				prot->p_topol->psi[r-1].num_fixed = 0;
				prot->p_topol->psi[r-1].fixed_atoms = NULL;
				prot->p_topol->psi[r-1].num_moved = 0;
				prot->p_topol->psi[r-1].moved_atoms = NULL;
				do_residue = 0;
			}
		}
		if (do_residue == 1){
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
	}
	
	free(atom_CA);
	free(atom_C);
	free(atom_O);

}

static void build_topology_individual_phi(protein_t *prot){
	char *atom_N, *atom_C, *atom_CA, *atom_O, *atom_H;
	int exists_atom_H;
	int i_af;
	int num_moved;
	int do_residue;

	atom_N = Malloc(char, 2);
	atom_CA = Malloc(char, 3);
	atom_C = Malloc(char, 2);	
	atom_O = Malloc(char, 2);
	atom_H = Malloc(char, 2);

	strcpy(atom_N, "N");
	strcpy(atom_CA, "CA");
	strcpy(atom_C, "C");
	strcpy(atom_O, "O");
	strcpy(atom_H, "H");

	//The first residue does not make rotation
	for (int r = 1; r <= prot->p_topol->numres; r++){		
		if (r > 1){
			do_residue = 1;
			if (r == prot->p_topol->numres){
				/* When residue is C-Termninal have is necessary check 
				* residue is ACE or NME. These residues are special.  
				*/
				if ( (strcmp(prot->p_atoms[0].resname, "ACE") == 0) ||
					 (strcmp(prot->p_atoms[0].resname, "NME") == 0) ){
					//Build Fixed and moved Atoms for ACE and NME
					prot->p_topol->phi[r-1].num_fixed = 0;
					prot->p_topol->phi[r-1].fixed_atoms = NULL;
					prot->p_topol->phi[r-1].num_moved = 0;
					prot->p_topol->phi[r-1].moved_atoms = NULL;
					do_residue = 0;
				}
			}
			if (do_residue == 1){
				exists_atom_H = 0;
				//Build Fixed Atoms
				num_moved = 2;
				i_af = -1;
				//num_hydrogen_backbone = get_number_hydrogen_backbone(prot, &r);
				if (atom_name_exists_in_resnum(prot->p_atoms, &r, 
					atom_H, &prot->p_topol->numatom) == btrue){				
					num_moved = num_moved + 1;
					exists_atom_H = 1;
				}
				prot->p_topol->phi[r-1].num_fixed = num_moved;
				prot->p_topol->phi[r-1].fixed_atoms = Malloc(int, prot->p_topol->phi[r-1].num_fixed);
				i_af++;
				prot->p_topol->phi[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
						&r, atom_N, &prot->p_topol->numatom);
				i_af++;
				prot->p_topol->phi[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
						&r, atom_CA, &prot->p_topol->numatom);
				if (exists_atom_H == 1){
					i_af++;
					prot->p_topol->phi[r-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
							&r, atom_H, &prot->p_topol->numatom);
				}

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
	free(atom_H);
}

/** Assigns the fixed and moved atoms at side chain based on residue name
* prot protein that will build the side chain topology
* res_num number of residue
* res_name means the name of residue
* chi represents what chi of residue 
*/
static void set_fixed_moved_atoms_side_chains_chi(protein_t *prot, 
	const int *res_num, const char *res_name, const int *chi){

	if ( strcmp(res_name, "SER") == 0){
		set_fixed_moved_atoms_side_chains_SER(prot, res_num, chi);
	}else if ( strcmp(res_name, "CYS") == 0){
		set_fixed_moved_atoms_side_chains_CYS(prot, res_num, chi);
	}else if ( strcmp(res_name, "THR") == 0){
		set_fixed_moved_atoms_side_chains_THR(prot, res_num, chi);
	}else if ( strcmp(res_name, "VAL") == 0){
		set_fixed_moved_atoms_side_chains_VAL(prot, res_num, chi);
	}else if ( strcmp(res_name, "ASP") == 0){
		set_fixed_moved_atoms_side_chains_ASP(prot, res_num, chi);
	}else if ( strcmp(res_name, "ASN") == 0){
		set_fixed_moved_atoms_side_chains_ASN(prot, res_num, chi);
	}else if ( strcmp(res_name, "ILE") == 0){
		set_fixed_moved_atoms_side_chains_ILE(prot, res_num, chi);
	}else if ( strcmp(res_name, "LEU") == 0){
		set_fixed_moved_atoms_side_chains_LEU(prot, res_num, chi);
	}else if ( strcmp(res_name, "PHE") == 0){
		set_fixed_moved_atoms_side_chains_PHE(prot, res_num, chi);
	}else if ( strcmp(res_name, "HIS") == 0){
		set_fixed_moved_atoms_side_chains_HIS(prot, res_num, chi);
	}else if ( strcmp(res_name, "TYR") == 0){
		set_fixed_moved_atoms_side_chains_TYR(prot, res_num, chi);
	}else if ( strcmp(res_name, "TRP") == 0){
		set_fixed_moved_atoms_side_chains_TRP(prot, res_num, chi);
	}else if ( strcmp(res_name, "MET") == 0){
		set_fixed_moved_atoms_side_chains_MET(prot, res_num, chi);
	}else if ( strcmp(res_name, "GLN") == 0){
		set_fixed_moved_atoms_side_chains_GLN(prot, res_num, chi);
	}else if ( strcmp(res_name, "GLU") == 0){
		set_fixed_moved_atoms_side_chains_GLU(prot, res_num, chi);
	}else if ( strcmp(res_name, "LYS") == 0){
		set_fixed_moved_atoms_side_chains_LYS(prot, res_num, chi);
	}else if ( strcmp(res_name, "ARG") == 0){
		set_fixed_moved_atoms_side_chains_ARG(prot, res_num, chi);
	}else{
		char msg[200];
		sprintf(msg," In set_fixed_moved_atoms_side_chains_chi function residue name %s was not found", res_name);
		fatal_error(msg);
	}
}

static void build_topology_individual_side_chains(protein_t *prot){
	char *res_name;
	res_name = Malloc(char,4);
	for (int r = 1; r <= prot->p_topol->numres; r++){
		get_res_name_from_res_num(res_name, &r, prot->p_atoms, &prot->p_topol->numatom);
		prot->p_topol->side_chains[r-1].num_chi = get_number_chi(res_name);
		if (prot->p_topol->side_chains[r-1].num_chi > 0){
			prot->p_topol->side_chains[r-1].atoms_chi = Malloc(top_residue_atom_info_t, 
				prot->p_topol->side_chains[r-1].num_chi);							
			for (int chi = 1; chi <= prot->p_topol->side_chains[r-1].num_chi; chi++){
				set_fixed_moved_atoms_side_chains_chi(prot, &r, res_name, &chi);				
			}
		}
	}
	free(res_name);
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



/** Builds the topology of protein
* prot is the protein building the topology
*/
_2PG_CARTESIAN_EXPORT
void build_topology_individual(protein_t *prot){
	build_topology_individual_atoms(prot);
	build_topology_individual_psi(prot);
	build_topology_individual_phi(prot);
	build_topology_individual_omega(prot);
	build_topology_individual_side_chains(prot);
}

/** Builds the topology of each protein (individual) which composes the population
* pop means the population
* pop_size size of population 
*/
_2PG_CARTESIAN_EXPORT
void build_topology_population(protein_t *pop, const int *pop_size){
	for (int i = 0; i < *pop_size; i++){
		build_topology_individual(&pop[i]);
	}
}

/** Returns the number of Hydrogens at backbone
*/
_2PG_CARTESIAN_EXPORT
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
_2PG_CARTESIAN_EXPORT
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
_2PG_CARTESIAN_EXPORT
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
_2PG_CARTESIAN_EXPORT
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
_2PG_CARTESIAN_EXPORT
boolean_t is_backbone_atom(const char *atomname){
	for (int i = 0; i < NUM_PROTEIN_BACKBONE_ATOMS ; i++){
		if ( strncmp(atomname, protein_backbone[i], 2) == 0){
			return btrue;
		}
	}
	return bfalse;
}

/** Returns the number of chi basead on residue name
* res_name is the residue name that wants to know the number
* of chi
*/
_2PG_CARTESIAN_EXPORT
int get_number_chi(const char *res_name){
		
	if( strcmp(res_name,"SER")==0 ){
		return 1;
	}else if( strcmp(res_name,"CYS")==0 ){
		return 1;	
	}else if( strcmp(res_name,"THR")==0 ){
		return 1;	
	}else if( strcmp(res_name,"VAL")==0 ){
		return 1;	
	}else if( strcmp(res_name,"ASP")==0 ){
		return 2;	
	}else if( strcmp(res_name,"ASN")==0 ){
		return 2;
	}else if( strcmp(res_name,"ILE")==0 ){
		return 2;
	}else if( strcmp(res_name,"LEU")==0 ){
		return 2;
	}else if( strcmp(res_name,"PRO")==0 ){
		return 0; //It is considered without chi, because it will be not moved.
	}else if( strcmp(res_name,"PHE")==0 ){
		return 2;
	}else if( strcmp(res_name,"HIS")==0 ){
		return 2;
	}else if( strcmp(res_name,"TYR")==0 ){
		return 2;
	}else if( strcmp(res_name,"TRP")==0 ){
		return 2;
	}else if( strcmp(res_name,"MET")==0 ){
		return 3;
	}else if( strcmp(res_name,"GLN")==0 ){
		return 3;
	}else if( strcmp(res_name,"GLU")==0 ){
		return 3;
	}else if( strcmp(res_name,"LYS")==0 ){
		return 4;
	}else if( strcmp(res_name,"ARG")==0 ){
		return 4;
	}else{
		return 0;
	}
}

/** Rename oxygen atoms in C-Terminal
* It is necessary because each force field adds two Oxygen atoms
* with diferent names. In Charmm27, the atoms are OT1 and OT2. In 
* Amber, the names are OC1 and OC2.

* atoms is the atoms
* res_num number of C-Terminal
* num_atom is the number of atoms
*/
_2PG_CARTESIAN_EXPORT
void rename_oxygen_c_terminal(pdb_atom_t *atoms,
		const int *res_num, const int *num_atom){

	char *atm_OT1, *atm_OT2, *atm_OC1, *atm_OC2, *atm_O;
	type_atoms t_OT1, t_OT2, t_OC1, t_OC2, t_O;
	int n_test;

	atm_OT1 = Malloc(char, 4);	
	atm_OT2 = Malloc(char, 4);
	atm_OC1 = Malloc(char, 4);
	atm_OC2 = Malloc(char, 4);
	atm_O = Malloc(char, 2);

	strcpy(atm_OT1, "OT1");
	strcpy(atm_OT2, "OT2");
	strcpy(atm_OC1, "OC1");
	strcpy(atm_OC2, "OC2");
	strcpy(atm_O, "O");

	t_OT1  = atmOT1;
	t_OT2  = atmOT2;
	t_OC1  = atmOC1; 
	t_OC2  = atmOC2; 
	t_O    = atmO; 

	// Is necessary look for two Oxygen atoms
	n_test = 1;	
	while (n_test <= 2){
		if (atom_name_exists_in_resnum(atoms, res_num, atm_OT1, num_atom) == btrue){
			rename_atom(atoms, atm_OT1, atm_O, res_num, &t_OT1, &t_O, num_atom);
		}else if (atom_name_exists_in_resnum(atoms, res_num, atm_OT2, num_atom) == btrue){
			rename_atom(atoms, atm_OT2, atm_O, res_num, &t_OT2, &t_O, num_atom);
		}else if (atom_name_exists_in_resnum(atoms, res_num, atm_OC1, num_atom) == btrue){
			rename_atom(atoms, atm_OC1, atm_O, res_num, &t_OC1, &t_O, num_atom);
		}else if (atom_name_exists_in_resnum(atoms, res_num, atm_OC2, num_atom) == btrue){
			rename_atom(atoms, atm_OC2, atm_O, res_num, &t_OC2, &t_O, num_atom);
		}
		n_test = n_test + 1;
	}	

	free(atm_OT1);
	free(atm_OT2);
	free(atm_OC1);
	free(atm_OC2);
	free(atm_O);

}
