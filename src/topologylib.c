#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "topologylib.h"
#include "messages.h"
#include "defines.h"
#include "pdbatom.h"

/** Returns the number of atom. If atomname does not exists
* a fatal_error is called. If you want to check atomname in 
* res_num use atom_name_exists_in_resnum function.
*
* atoms of molecule
* res_num number of residue in wich want to check atom 
* atomname name of atom that want to check if exists in res_num
* num_atom number of atoms of molecule
*/
int get_atom_index_by_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname, const int *num_atom){
	const pdb_atom_t * aux;

	aux = search_pdb_atom_from_resnum_atom_name(atoms, res_num, atomname, num_atom);
	if (aux == NULL){
 		char msg[300];
 	    sprintf(msg,"Atom %s not found at get_atom_index_by_resnum_atom_name \n", atomname);
  		fatal_error(msg);		
	}
	return aux->atmnumber;
}

/** checks if atom name exists in residue. 
* Returns true when exists. Otherwise, returns false.
*
* atoms of molecule
* res_num number of residue in wich want to check atom 
* atomname name of atom that want to check if exists in res_num
* num_atom number of atoms of molecule
*/
boolean_t atom_name_exists_in_resnum(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname, const int *num_atom){
	const pdb_atom_t * aux;

	aux = search_pdb_atom_from_resnum_atom_name(atoms, res_num, atomname, num_atom);
	if (aux != NULL){
		return btrue;
	}
	return bfalse;
}

/** Assigens the fixed and moved atoms for SER
*/
void set_fixed_moved_atoms_side_chains_SER(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_OG, *atom_HG1, *atom_HB1, *atom_HB2;
	int i_af;
	int num_moved;
	int exists_HG1, exists_HB1, exists_HB2;

	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_OG = Malloc(char, 3);
	atom_HG1 = Malloc(char, 4);
	atom_HB1 = Malloc(char, 4);
	atom_HB2 = Malloc(char, 4);

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_OG, "OG");
	strcpy(atom_HG1, "HG1");	
	strcpy(atom_HB1, "HB1");
	strcpy(atom_HB2, "HB2");

	exists_HG1 = 0;
	exists_HB1 = 0;
	exists_HB2 = 0;

	//Building Fixed Atoms
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CA, &prot->p_topol->numatom);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CB, &prot->p_topol->numatom);
	//Building Moved Atoms
	num_moved = 1;
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_HB1 = 1;
	}		
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_HB2 = 1;
	}			
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_HG1 = 1;
	}	
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_OG, &prot->p_topol->numatom);
	if (exists_HB1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom);
	}	
	if (exists_HB2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HB2, &prot->p_topol->numatom);
	}		
	if (exists_HG1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HG1, &prot->p_topol->numatom);
	}
	free(atom_CA);
	free(atom_CB);
	free(atom_OG);
	free(atom_HG1);
	free(atom_HB1);
	free(atom_HB2);	
}

/** Assigens the fixed and moved atoms for CYS
*/
void set_fixed_moved_atoms_side_chains_CYS(protein_t *prot, 
	const int *res_num, const int *chi){	
	char *atom_CA, *atom_CB, *atom_SG, *atom_HG1, *atom_HB1, *atom_HB2;
	int exists_atom_HG1, exists_HB1, exists_HB2;
	int num_moved;
	int i_af;

	exists_atom_HG1 = 0;
	exists_HB1 = 0;
	exists_HB2 = 0;
	
	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_SG = Malloc(char, 3);
	atom_HG1 = Malloc(char, 4);
	atom_HB1 = Malloc(char, 4);
	atom_HB2 = Malloc(char, 4);

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_SG, "SG");
	strcpy(atom_HG1, "HG1");
	strcpy(atom_HB1, "HB1");
	strcpy(atom_HB2, "HB2");

	//Building Fixed Atoms
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CA, &prot->p_topol->numatom);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CB, &prot->p_topol->numatom);
	//Building Moved Atoms
	num_moved = 1;
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_HB1 = 1;
	}		
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_HB2 = 1;
	}	
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_HG1 = 1;
	}	
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_SG, &prot->p_topol->numatom);
	if (exists_HB1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom);
	}	
	if (exists_HB2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HB2, &prot->p_topol->numatom);
	}			
	if (exists_atom_HG1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HG1, &prot->p_topol->numatom);

	}
	free(atom_CA);
	free(atom_CB);
	free(atom_SG);
	free(atom_HG1);
	free(atom_HB1);
	free(atom_HB2);	

}

/** Assigens the fixed and moved atoms for THR
*/
void set_fixed_moved_atoms_side_chains_THR(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_OG1, *atom_CG2, 
	*atom_HG1, *atom_1HG2, *atom_2HG2, *atom_3HG2, *atom_HB;
	int i_af;
	int num_moved;
	int exists_atom_HG1, exists_atom_1HG2, 
	exists_atom_2HG2, exists_atom_3HG2, exists_atom_HB;
	
	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_OG1 = Malloc(char, 4);
	atom_CG2 = Malloc(char, 4);
	atom_HG1  = Malloc(char, 4); 
	atom_1HG2 = Malloc(char, 5); 
	atom_2HG2 = Malloc(char, 5);
	atom_3HG2  = Malloc(char, 5);
	atom_HB = Malloc(char, 3);	

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_OG1, "OG1");
	strcpy(atom_CG2, "CG2");
	strcpy(atom_HG1, "HG1");
	strcpy(atom_1HG2, "1HG2");
	strcpy(atom_2HG2, "2HG2");
	strcpy(atom_3HG2, "3HG2");
	strcpy(atom_HB, "HB");	
	
	exists_atom_HG1 = 0;
	exists_atom_1HG2 = 0;
	exists_atom_2HG2 = 0; 
	exists_atom_3HG2 = 0;
	exists_atom_HB = 0;

	//Building Fixed Atoms
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CA, &prot->p_topol->numatom);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CB, &prot->p_topol->numatom);
	//Building Moved Atoms
	num_moved = 2; 
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_HB = 1;
	}	
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_HG1 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HG2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_1HG2 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HG2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_2HG2 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_3HG2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_3HG2 = 1;
	}
	i_af = -1;	
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_OG1, &prot->p_topol->numatom);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CG2, &prot->p_topol->numatom);
	if (exists_atom_HB == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HB, &prot->p_topol->numatom);
	}	
	if (exists_atom_HG1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom);
	}
	if (exists_atom_1HG2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_1HG2, &prot->p_topol->numatom);
	}
	if (exists_atom_2HG2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_2HG2, &prot->p_topol->numatom);
	}
	if (exists_atom_3HG2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_3HG2, &prot->p_topol->numatom);
	}

	free(atom_CA);
	free(atom_CB);
	free(atom_OG1);
	free(atom_CG2);
	free(atom_HG1);
	free(atom_1HG2);
	free(atom_2HG2);
	free(atom_3HG2);
	free(atom_HB);
}

/** Assigens the fixed and moved atoms for VAL
*/
void set_fixed_moved_atoms_side_chains_VAL(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG1, *atom_CG2, 
	*atom_1HG1, *atom_2HG1, *atom_3HG1, *atom_1HG2, 
	*atom_2HG2, *atom_3HG2, *atom_HB;
	int i_af;
	int exists_atom_1HG1, exists_atom_2HG1, 
	exists_atom_3HG1, exists_atom_1HG2, exists_atom_2HG2, 
	exists_atom_3HG2, exists_atom_HB;
	int num_moved;
	
	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_CG1 = Malloc(char, 4);
	atom_CG2 = Malloc(char, 4);	
	atom_1HG1 = Malloc(char, 5); 
	atom_2HG1 = Malloc(char, 5); 
	atom_3HG1 = Malloc(char, 5);
	atom_1HG2 = Malloc(char, 5);
	atom_2HG2 = Malloc(char, 5); 
	atom_3HG2 = Malloc(char, 5);
	atom_HB = Malloc(char, 3);

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_CG1, "CG1");
	strcpy(atom_CG2, "CG2");
	strcpy(atom_1HG1, "1HG1");
	strcpy(atom_2HG1, "2HG1");
	strcpy(atom_3HG1, "3HG1");
	strcpy(atom_1HG2, "1HG2");
	strcpy(atom_2HG2, "2HG2");
	strcpy(atom_3HG2, "3HG2");
	strcpy(atom_HB, "HB");

	exists_atom_1HG1  = 0; 
	exists_atom_2HG1  = 0;
	exists_atom_3HG1  = 0;
	exists_atom_1HG2  = 0; 
	exists_atom_2HG2  = 0; 
	exists_atom_3HG2  = 0;
	exists_atom_HB    = 0;

	//Building Fixed Atoms
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CA, &prot->p_topol->numatom);
	i_af++;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CB, &prot->p_topol->numatom);
	//Building Moved Atoms
	num_moved = 2; 
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_HB = 1;
	}	
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HG1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_1HG1 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HG1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_2HG1 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_3HG1, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_3HG1 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HG2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_1HG2 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HG2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_2HG2 = 1;
	}
	if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_3HG2, &prot->p_topol->numatom) == btrue){
		num_moved = num_moved + 1;
		exists_atom_3HG2 = 1;
	}
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CG1, &prot->p_topol->numatom);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CG2, &prot->p_topol->numatom);
	if (exists_atom_HB == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_HB, &prot->p_topol->numatom);		
	}
	if (exists_atom_1HG1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_1HG1, &prot->p_topol->numatom);		
	}
	if (exists_atom_2HG1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_2HG1, &prot->p_topol->numatom);		
	}
	if (exists_atom_3HG1 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_3HG1, &prot->p_topol->numatom);		
	}
	if (exists_atom_1HG2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_1HG2, &prot->p_topol->numatom);		
	}
	if (exists_atom_2HG2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_2HG2, &prot->p_topol->numatom);		
	}
	if (exists_atom_3HG2 == 1){
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_3HG2, &prot->p_topol->numatom);		
	}
	free(atom_CA);
	free(atom_CB);
	free(atom_CG1);
	free(atom_CG2);
	free(atom_HB);
	free(atom_1HG1);
	free(atom_2HG1);
	free(atom_3HG1);
	free(atom_1HG2);
	free(atom_2HG2);
	free(atom_3HG2);

}

/** Assigens the fixed and moved atoms for ASP
*/
void set_fixed_moved_atoms_side_chains_ASP(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_OD1, *atom_OD2,
	*atom_HD2, *atom_HB1, *atom_HB2;
	int i_af;	
	int exists_atom_HD2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	
	exists_atom_HD2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_OD2 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_OD2, "OD2");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}			
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_OD2);		
		free(atom_HD2);
		free(atom_HB1);
		free(atom_HB2);			

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_OD2 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);


		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_OD2, "OD2");
		strcpy(atom_HD2, "HD2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 2;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD2, &prot->p_topol->numatom);
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		
		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_OD2);
		free(atom_HD2);

	}

}

/** Assigens the fixed and moved atoms for ASN
*/
void set_fixed_moved_atoms_side_chains_ASN(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_OD1, *atom_ND2, 
	*atom_1HD2, *atom_2HD2, *atom_HB1, *atom_HB2;
	int exists_atom_1HD2, exists_atom_2HD2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_1HD2 = 0;
	exists_atom_2HD2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_ND2 = Malloc(char, 4);
		atom_1HD2 = Malloc(char, 5);
		atom_2HD2 = Malloc(char, 5);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_ND2, "ND2");
		strcpy(atom_1HD2, "1HD2");
		strcpy(atom_2HD2, "2HD2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HD2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_ND2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_1HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_2HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HD2, &prot->p_topol->numatom);
		}
		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_ND2);
		free(atom_1HD2);
		free(atom_2HD2);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_ND2 = Malloc(char, 4);
		atom_1HD2 = Malloc(char, 5);
		atom_2HD2 = Malloc(char, 5);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_ND2, "ND2");
		strcpy(atom_1HD2, "1HD2");
		strcpy(atom_2HD2, "2HD2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 2;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HD2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_ND2, &prot->p_topol->numatom);
		if (exists_atom_1HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_2HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HD2, &prot->p_topol->numatom);
		}		

		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_ND2);
		free(atom_1HD2);
		free(atom_2HD2);

	}

}

/** Assigens the fixed and moved atoms for ILE
*/
void set_fixed_moved_atoms_side_chains_ILE(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_HB, *atom_CG2, *atom_1HG2, 
	*atom_2HG2, *atom_3HG2, *atom_CG1, *atom_1HG1, *atom_2HG1, 
	*atom_CD, *atom_HD1, *atom_HD2, *atom_HD3;
	int i_af;
	int exists_atom_1HG1, exists_atom_2HG1, 
	exists_atom_1HG2, exists_atom_2HG2, 
	exists_atom_3HG2, exists_atom_HB, exists_atom_HD1,
	exists_atom_HD2, exists_atom_HD3;
	int num_moved;

	exists_atom_1HG1  = 0; 
	exists_atom_2HG1  = 0;
	exists_atom_1HG2  = 0; 
	exists_atom_2HG2  = 0; 
	exists_atom_3HG2  = 0;
	exists_atom_HB    = 0;
	exists_atom_HD1   = 0;
	exists_atom_HD2   = 0;
	exists_atom_HD3	  = 0;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG1 = Malloc(char, 4);
		atom_CG2 = Malloc(char, 4);	
		atom_1HG1 = Malloc(char, 5); 
		atom_2HG1 = Malloc(char, 5); 		
		atom_1HG2 = Malloc(char, 5);
		atom_2HG2 = Malloc(char, 5); 
		atom_3HG2 = Malloc(char, 5);
		atom_HB = Malloc(char, 3);
		atom_CD =  Malloc(char, 3);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HD3 = Malloc(char, 4);		

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG1, "CG1");
		strcpy(atom_CD, "CD");		
		strcpy(atom_CG2, "CG2");
		strcpy(atom_HB, "HB");

		strcpy(atom_1HG1, "1HG1");
		strcpy(atom_2HG1, "2HG1");		
		strcpy(atom_1HG2, "1HG2");
		strcpy(atom_2HG2, "2HG2");
		strcpy(atom_3HG2, "3HG2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HD3, "HD3");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3; 
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_1HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HG1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_2HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HG1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD3 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_1HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HG2 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_2HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HG2 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_3HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_3HG2 = 1;
		}

		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG1, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG2, &prot->p_topol->numatom);

		if (exists_atom_HB == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB, &prot->p_topol->numatom);		
		}
		if (exists_atom_1HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HG1, &prot->p_topol->numatom);		
		}
		if (exists_atom_2HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HG1, &prot->p_topol->numatom);		
		}
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);		
		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);		
		}
		if (exists_atom_HD3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD3, &prot->p_topol->numatom);		
		}
		if (exists_atom_1HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HG2, &prot->p_topol->numatom);		
		}
		if (exists_atom_2HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HG2, &prot->p_topol->numatom);		
		}
		if (exists_atom_3HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_3HG2, &prot->p_topol->numatom);		
		}
		
		free(atom_CA);
		free(atom_CB);
		free(atom_CG1);
		free(atom_CG2);
		free(atom_HB);
		free(atom_1HG1);
		free(atom_2HG1);
		free(atom_1HG2);
		free(atom_2HG2);
		free(atom_3HG2);
		free(atom_CD);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HD3);

	}else if (*chi == 2){		
		
		atom_CB = Malloc(char, 3);
		atom_CD = Malloc(char, 3);			
		atom_CG1 = Malloc(char, 4);		
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HD3 = Malloc(char, 4);		
		atom_1HG1 = Malloc(char, 5); 
		atom_2HG1 = Malloc(char, 5); 		
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CD, "CD");
		strcpy(atom_CG1, "CG1");		
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HD3, "HD3");
		strcpy(atom_1HG1, "1HG1");
		strcpy(atom_2HG1, "2HG1");		

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG1, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 1; 
		
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD3 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_1HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HG1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_2HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HG1 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);

		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);		
		}		
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);		
		}		
		if (exists_atom_HD3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD3, &prot->p_topol->numatom);		
		}
		if (exists_atom_1HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HG1, &prot->p_topol->numatom);		
		}
		if (exists_atom_2HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HG1, &prot->p_topol->numatom);		
		}

		free(atom_CB);
		free(atom_CD);
		free(atom_CG1);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HD3);
		free(atom_1HG1);
		free(atom_2HG1);
		
	}
}

/** Assigens the fixed and moved atoms for LEU
*/
void set_fixed_moved_atoms_side_chains_LEU(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2, 
	*atom_1HD1, *atom_2HD1, *atom_3HD1, *atom_1HD2, *atom_2HD2, 
	*atom_3HD2, *atom_HB1, *atom_HB2, *atom_HG;
	int exists_atom_1HD1, exists_atom_2HD1, exists_atom_3HD1, 
	exists_atom_1HD2, exists_atom_2HD2, exists_atom_3HD2, 
	exists_atom_HB1, exists_atom_HB2, exists_atom_HG;
	int num_moved;
	int i_af;

	exists_atom_1HD1 = 0;
	exists_atom_2HD1 = 0;
	exists_atom_3HD1 = 0; 
	exists_atom_1HD2 = 0; 
	exists_atom_2HD2 = 0;
	exists_atom_3HD2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;
	exists_atom_HG = 0;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_1HD1 = Malloc(char, 5);
		atom_2HD1 = Malloc(char, 5);
		atom_3HD1 = Malloc(char, 5);
		atom_1HD2 = Malloc(char, 5);
		atom_2HD2 = Malloc(char, 5);
		atom_3HD2 = Malloc(char, 5);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);
		atom_HG = Malloc(char, 3);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1, "CD1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_1HD1, "1HD1");
		strcpy(atom_2HD1, "2HD1");
		strcpy(atom_3HD1, "3HD1");
		strcpy(atom_1HD2, "1HD2");
		strcpy(atom_2HD2, "2HD2");
		strcpy(atom_3HD2, "3HD2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");
		strcpy(atom_HG, "HG");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}			
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG = 1;
		}					
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_3HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_3HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_3HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_3HD2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}
		if (exists_atom_HG == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG, &prot->p_topol->numatom);
		}		
		if (exists_atom_1HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_2HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_3HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_3HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_1HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_2HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_3HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_3HD2, &prot->p_topol->numatom);
		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_1HD1);
		free(atom_2HD1);
		free(atom_3HD1);
		free(atom_1HD2);
		free(atom_2HD2);
		free(atom_3HD2);
		free(atom_HB1);
		free(atom_HB2);	
		free(atom_HG);	
	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_1HD1 = Malloc(char, 5);
		atom_2HD1 = Malloc(char, 5);
		atom_3HD1 = Malloc(char, 5);
		atom_1HD2 = Malloc(char, 5);
		atom_2HD2 = Malloc(char, 5);
		atom_3HD2 = Malloc(char, 5);
		atom_HG = Malloc(char, 3);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1, "CD1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_1HD1, "1HD1");
		strcpy(atom_2HD1, "2HD1");
		strcpy(atom_3HD1, "3HD1");
		strcpy(atom_1HD2, "1HD2");
		strcpy(atom_2HD2, "2HD2");
		strcpy(atom_3HD2, "3HD2");
		strcpy(atom_HG, "HG");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 2;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_3HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_3HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_3HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_3HD2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		if (exists_atom_HG == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG, &prot->p_topol->numatom);
		}
		if (exists_atom_1HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_2HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_3HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_3HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_1HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_2HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_3HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_3HD2, &prot->p_topol->numatom);
		}
		
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_1HD1);
		free(atom_2HD1);
		free(atom_3HD1);
		free(atom_1HD2);
		free(atom_2HD2);
		free(atom_3HD2);
		free(atom_HG);

	}

}



/** Assigens the fixed and moved atoms for PHE
*/
void set_fixed_moved_atoms_side_chains_PHE(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2, 
	*atom_CE1, *atom_CE2, *atom_CZ, *atom_HD1, *atom_HE1, *atom_HZ, 
	*atom_HD2, *atom_HE2, *atom_HB1, *atom_HB2;
	int exists_atom_HD1, exists_atom_HE1, exists_atom_HZ, exists_atom_HD2, 
	exists_atom_HE2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_HD1 = 0; 
	exists_atom_HE1 = 0; 
	exists_atom_HZ = 0; 
	exists_atom_HD2 = 0; 
	exists_atom_HE2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ = Malloc(char, 3);
		atom_HD1 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HZ = Malloc(char, 3);
		atom_HD2 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ, "CZ" );
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HZ, "HZ");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 6;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);
		free(atom_HD1); 
		free(atom_HE1); 
		free(atom_HZ); 
		free(atom_HD2); 
		free(atom_HE2);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ = Malloc(char, 3);
		atom_HD1 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HZ = Malloc(char, 3);
		atom_HD2 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ,"CZ" );
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HZ, "HZ");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE2, "HE2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved =  5;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}

		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);
		free(atom_HD1); 
		free(atom_HE1); 
		free(atom_HZ); 
		free(atom_HD2); 
		free(atom_HE2);

		
	}

}


/** Assigens the fixed and moved atoms for HIS
*/
void set_fixed_moved_atoms_side_chains_HIS(protein_t *prot, 
	const int *res_num, const int *chi){
	char *atom_CA, *atom_CB, *atom_CG, *atom_ND1, *atom_CD2, *atom_CE1, *atom_NE2, 
	*atom_HD1, *atom_HE1, *atom_HD2, *atom_HE2, *atom_HB1, *atom_HB2;
	int exists_atom_HD1, exists_atom_HE1, exists_atom_HD2, 
	exists_atom_HE2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_HD1 = 0; 
	exists_atom_HD2 = 0; 	
	exists_atom_HE1 = 0; 	
	exists_atom_HE2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_ND1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);	
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_ND1, "ND1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_CE1, "CE1");
		strcpy(atom_NE2, "NE2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved =	5;	
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_ND1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);
		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);
		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_ND1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_NE2);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HE1);
		free(atom_HE2);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_ND1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);	
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_ND1, "ND1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_CE1, "CE1");
		strcpy(atom_NE2, "NE2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 4;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_ND1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE2, &prot->p_topol->numatom);
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);
		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);
		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);
		}

		free(atom_CB);
		free(atom_CG);
		free(atom_ND1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_NE2);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HE1);
		free(atom_HE2);
		
	}

}


/** Assigens the fixed and moved atoms for TYR
*/
void set_fixed_moved_atoms_side_chains_TYR(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2, 
	*atom_CE1, *atom_CE2, *atom_CZ, *atom_OH, *atom_HD1, *atom_HE1, 
	*atom_HH, *atom_HD2, *atom_HE2, *atom_HB1, *atom_HB2;
	int exists_atom_HD1, exists_atom_HE1, exists_atom_HH, 
	exists_atom_HD2, exists_atom_HE2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_HD1 = 0; 
	exists_atom_HE1 = 0; 
	exists_atom_HH = 0;
	exists_atom_HD2 = 0; 
	exists_atom_HE2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ = Malloc(char, 3);
		atom_OH = Malloc(char, 3);
		atom_HD1 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HH = Malloc(char, 3);
		atom_HD2 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ,"CZ" );
		strcpy(atom_OH,"OH" );
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HH, "HH");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 7;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HH, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HH = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		

		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OH, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HH == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HH, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);
		free(atom_OH);
		free(atom_HD1);
		free(atom_HE1);
		free(atom_HH);
		free(atom_HD2);
		free(atom_HE2);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ = Malloc(char, 3);
		atom_OH = Malloc(char, 3);
		atom_HD1 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HH = Malloc(char, 3);
		atom_HD2 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ,"CZ" );
		strcpy(atom_OH,"OH" );
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HH, "HH");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE2, "HE2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 6;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HH, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HH = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}				
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OH, &prot->p_topol->numatom);
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HH == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HH, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}


		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);
		free(atom_OH);
		free(atom_HD1);
		free(atom_HE1);
		free(atom_HH);
		free(atom_HD2);
		free(atom_HE2);

		
	}

}


/** Assigens the fixed and moved atoms for TRP
*/
void set_fixed_moved_atoms_side_chains_TRP(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2, *atom_NE1, *atom_CE3, 
	*atom_CE2, *atom_CZ3, *atom_CZ2, *atom_CH2, 
	*atom_HD1, *atom_HE1, *atom_HE3, *atom_HZ3, *atom_HZ2, 
	*atom_HH2, *atom_HB1, *atom_HB2;
	int exists_atom_HD1, exists_atom_HE1, exists_atom_HE3, exists_atom_HZ3, 
	exists_atom_HZ2, exists_atom_HH2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_HD1 = 0; 
	exists_atom_HE1 = 0; 
	exists_atom_HE3 = 0; 
	exists_atom_HZ3 = 0; 
	exists_atom_HZ2 = 0; 
	exists_atom_HH2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_NE1 = Malloc(char, 4);
		atom_CE3 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ3 = Malloc(char, 4);
		atom_CZ2 = Malloc(char, 4);
		atom_CH2 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE3 = Malloc(char, 4);
		atom_HZ3 = Malloc(char, 4);
		atom_HZ2 = Malloc(char, 4);
		atom_HH2 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1, "CD1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_NE1, "NE1");
		strcpy(atom_CE3, "CE3");
		strcpy(atom_CE2, "CE2");
		strcpy(atom_CZ3, "CZ3");
		strcpy(atom_CZ2, "CZ2");
		strcpy(atom_CH2, "CH2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE3, "HE3");
		strcpy(atom_HZ3, "HZ3");
		strcpy(atom_HZ2, "HZ2");
		strcpy(atom_HH2, "HH2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 9;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}						
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE3 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ3 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HH2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE3, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ3, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CH2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);
		}		
		if (exists_atom_HE3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE3, &prot->p_topol->numatom);
		}
		if (exists_atom_HZ3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ3, &prot->p_topol->numatom);
		}
		if (exists_atom_HZ2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ2, &prot->p_topol->numatom);
		}
		if (exists_atom_HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HH2, &prot->p_topol->numatom);
		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_NE1);
		free(atom_CE3);
		free(atom_CE2);
		free(atom_CZ3);
		free(atom_CZ2);
		free(atom_CH2);
		free(atom_HD1);
		free(atom_HE1);
		free(atom_HE3);
		free(atom_HZ3);
		free(atom_HZ2);
		free(atom_HH2);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_NE1 = Malloc(char, 4);
		atom_CE3 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ3 = Malloc(char, 4);
		atom_CZ2 = Malloc(char, 4);
		atom_CH2 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE3 = Malloc(char, 4);
		atom_HZ3 = Malloc(char, 4);
		atom_HZ2 = Malloc(char, 4);
		atom_HH2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1, "CD1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_NE1, "NE1");
		strcpy(atom_CE3, "CE3");
		strcpy(atom_CE2, "CE2");
		strcpy(atom_CZ3, "CZ3");
		strcpy(atom_CZ2, "CZ2");
		strcpy(atom_CH2, "CH2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE3, "HE3");
		strcpy(atom_HZ3, "HZ3");
		strcpy(atom_HZ2, "HZ2");
		strcpy(atom_HH2, "HH2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 8;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}						
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE3 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ3 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HH2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE3, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ3, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ2, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CH2, &prot->p_topol->numatom);
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);
		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);
		}		
		if (exists_atom_HE3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE3, &prot->p_topol->numatom);
		}
		if (exists_atom_HZ3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ3, &prot->p_topol->numatom);
		}
		if (exists_atom_HZ2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ2, &prot->p_topol->numatom);
		}
		if (exists_atom_HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HH2, &prot->p_topol->numatom);
		}

		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_NE1);
		free(atom_CE3);
		free(atom_CE2);
		free(atom_CZ3);
		free(atom_CZ2);
		free(atom_CH2);
		free(atom_HD1);
		free(atom_HE1);
		free(atom_HE3);
		free(atom_HZ3);
		free(atom_HZ2);
		free(atom_HH2);

	}

}

/** Assigens the fixed and moved atoms for MET
*/
void set_fixed_moved_atoms_side_chains_MET(protein_t *prot, 
	const int *res_num, const int *chi){
	char *atom_CA, *atom_CB, *atom_CG, *atom_SD, *atom_CE, 
	*atom_HG1, *atom_HG2, *atom_HE1, *atom_HE2, 
	*atom_HE3, *atom_HB1, *atom_HB2;
	int exists_atom_HG1, exists_atom_HG2, exists_atom_HE1, 
	exists_atom_HE2, exists_atom_HE3, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_HG1 = 0;
	exists_atom_HG2 = 0; 
	exists_atom_HE1 = 0;
	exists_atom_HE2 = 0; 
	exists_atom_HE3 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_SD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HE3 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_SD, "SD");
		strcpy(atom_CE, "CE");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HE3, "HE3");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE3 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_SD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);

		}
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE3, &prot->p_topol->numatom);

		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_SD);
		free(atom_CE);
		free(atom_HG1);
		free(atom_HG2);
		free(atom_HE1);
		free(atom_HE2);
		free(atom_HE3);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_SD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HE3 = Malloc(char, 4);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_SD, "SD");
		strcpy(atom_CE, "CE");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HE3, "HE3");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 2;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE3 = 1;
		}		

		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_SD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE3, &prot->p_topol->numatom);

		}		
	
		free(atom_CB);
		free(atom_CG);
		free(atom_SD);
		free(atom_CE);
		free(atom_HE1);
		free(atom_HE2);
		free(atom_HE3);
		free(atom_HG1);
		free(atom_HG2);			
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_SD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HE3 = Malloc(char, 4);
		
		strcpy(atom_CG, "CG");
		strcpy(atom_SD, "SD");
		strcpy(atom_CE, "CE");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HE3, "HE3");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_SD, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 1;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE3 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE3, &prot->p_topol->numatom);

		}		
			
		free(atom_CG);
		free(atom_SD);
		free(atom_CE);	
		free(atom_HE1);
		free(atom_HE2);
		free(atom_HE3);

	}

}

/** Assigens the fixed and moved atoms for GLN
*/
void set_fixed_moved_atoms_side_chains_GLN(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_OE1, *atom_NE2, 
	*atom_HG1, *atom_HG2, *atom_1HE2, *atom_2HE2, *atom_HB1, *atom_HB2;
	int exists_atom_HG1, exists_atom_HG2, exists_atom_1HE2, 
	exists_atom_2HE2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_HG1 = 0; 
	exists_atom_HG2 = 0; 
	exists_atom_1HE2 = 0; 
	exists_atom_2HE2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);
		atom_1HE2 = Malloc(char, 5);
		atom_2HE2 = Malloc(char, 5);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_NE2, "NE2");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");
		strcpy(atom_1HE2, "1HE2");
		strcpy(atom_2HE2, "2HE2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 4;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HE2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HE2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);

		}
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);

		}
		if (exists_atom_1HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_2HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HE2, &prot->p_topol->numatom);

		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_NE2);
		free(atom_HG1);
		free(atom_HG2);
		free(atom_1HE2);
		free(atom_2HE2);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
		atom_1HE2 = Malloc(char, 5);
		atom_2HE2 = Malloc(char, 5);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_NE2, "NE2");
		strcpy(atom_1HE2, "1HE2");
		strcpy(atom_2HE2, "2HE2");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;		
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HE2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HE2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE2, &prot->p_topol->numatom);
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);
		}
		if (exists_atom_1HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_2HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HE2, &prot->p_topol->numatom);

		}

		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_NE2);
		free(atom_1HE2);
		free(atom_2HE2);
		free(atom_HG1);
		free(atom_HG2);	
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
		atom_1HE2 = Malloc(char, 5);
		atom_2HE2 = Malloc(char, 5);
				
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_NE2, "NE2");
		strcpy(atom_1HE2, "1HE2");
		strcpy(atom_2HE2, "2HE2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 2;		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HE2 = 1;
		}				
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HE2 = 1;
		}				
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE2, &prot->p_topol->numatom);
		if (exists_atom_1HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_2HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HE2, &prot->p_topol->numatom);

		}

		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_NE2);
		free(atom_1HE2);
		free(atom_2HE2);

	}

}


/** Assigens the fixed and moved atoms for GLU
*/
void set_fixed_moved_atoms_side_chains_GLU(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_OE1, *atom_OE2, 
	*atom_HG1, *atom_HG2, *atom_HE2, *atom_HB1, *atom_HB2;
	int exists_atom_HG1, exists_atom_HG2, exists_atom_HE2, exists_atom_HB1, 
	exists_atom_HB2;
	int num_moved;
	int i_af;

	exists_atom_HG1 = 0; 
	exists_atom_HG2 = 0; 
	exists_atom_HE2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_OE2 = Malloc(char, 4);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_OE2, "OE2");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 4;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);
		}
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);
		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);
		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_OE2);
		free(atom_HG1);
		free(atom_HG2);
		free(atom_HE2);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_OE2 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_OE2, "OE2");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE2, &prot->p_topol->numatom);
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_OE2);
		free(atom_HE2);
		free(atom_HG1);
		free(atom_HG2);	
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_OE2 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
				
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_OE2, "OE2");
		strcpy(atom_HE2, "HE2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 2;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE2, &prot->p_topol->numatom);
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}

		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_OE2);
		free(atom_HE2);

	}

}


/** Assigens the fixed and moved atoms for LYS
*/
void set_fixed_moved_atoms_side_chains_LYS(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_CE, *atom_NZ, 
	*atom_HG1, *atom_HG2, *atom_HD1, *atom_HD2, *atom_HE1, *atom_HE2, 
	*atom_HZ1, *atom_HZ2, *atom_HZ3, *atom_HB1, *atom_HB2;
	int exists_atom_HG1, exists_atom_HG2, exists_atom_HD1, exists_atom_HD2, 
	exists_atom_HE1, exists_atom_HE2, exists_atom_HZ1, exists_atom_HZ2, 
	exists_atom_HZ3, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;
	
	exists_atom_HG1 = 0; 
	exists_atom_HG2 = 0; 
	exists_atom_HD1 = 0; 
	exists_atom_HD2 = 0;
	exists_atom_HE1 = 0; 
	exists_atom_HE2 = 0; 
	exists_atom_HZ1 = 0; 
	exists_atom_HZ2 = 0;
	exists_atom_HZ3 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char, 3);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HZ1 = Malloc(char, 4);
		atom_HZ2 = Malloc(char, 4);
		atom_HZ3 = Malloc(char, 4);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);
		
		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HZ1, "HZ1");
		strcpy(atom_HZ2, "HZ2");
		strcpy(atom_HZ3, "HZ3");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 4;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ3 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NZ, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);

		}
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);

		}
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ1, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ2, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ3, &prot->p_topol->numatom);

		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);
		free(atom_HG1);
		free(atom_HG2);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HE1);
		free(atom_HE2);
		free(atom_HZ1);
		free(atom_HZ2);
		free(atom_HZ3);
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char, 3);		
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HZ1 = Malloc(char, 4);
		atom_HZ2 = Malloc(char, 4);
		atom_HZ3 = Malloc(char, 4);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HZ1, "HZ1");
		strcpy(atom_HZ2, "HZ2");
		strcpy(atom_HZ3, "HZ3");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ3 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NZ, &prot->p_topol->numatom);
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);
		}
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ1, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ2, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ3, &prot->p_topol->numatom);

		}

		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HE1);
		free(atom_HE2);
		free(atom_HZ1);
		free(atom_HZ2);
		free(atom_HZ3);
		free(atom_HG1);
		free(atom_HG2);	
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char,3);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
		atom_HZ1 = Malloc(char, 4);
		atom_HZ2 = Malloc(char, 4);
		atom_HZ3 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
				
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");
		strcpy(atom_HZ1, "HZ1");
		strcpy(atom_HZ2, "HZ2");
		strcpy(atom_HZ3, "HZ3");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 2;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ3 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NZ, &prot->p_topol->numatom);
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ1, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ2, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ3, &prot->p_topol->numatom);

		}

		free(atom_CG);
		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);
		free(atom_HE1);
		free(atom_HE2);
		free(atom_HZ1);
		free(atom_HZ2);
		free(atom_HZ3);
		free(atom_HD1);
		free(atom_HD2);	

	}else if (*chi == 4){
				
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char,3);
		atom_HZ1 = Malloc(char, 4);
		atom_HZ2 = Malloc(char, 4);
		atom_HZ3 = Malloc(char, 4);
		atom_HE1 = Malloc(char, 4);
		atom_HE2 = Malloc(char, 4);
						
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");
		strcpy(atom_HZ1, "HZ1");
		strcpy(atom_HZ2, "HZ2");
		strcpy(atom_HZ3, "HZ3");
		strcpy(atom_HE1, "HE1");
		strcpy(atom_HE2, "HE2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 1;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HZ3, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HZ3 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE2 = 1;
		}
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NZ, &prot->p_topol->numatom);
		if (exists_atom_HZ1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ1, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ2, &prot->p_topol->numatom);

		}
		if (exists_atom_HZ3 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HZ3, &prot->p_topol->numatom);

		}
		if (exists_atom_HE1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE1, &prot->p_topol->numatom);

		}		
		if (exists_atom_HE2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE2, &prot->p_topol->numatom);

		}				
		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);
		free(atom_HZ1);
		free(atom_HZ2);
		free(atom_HZ3);
		free(atom_HE1);
		free(atom_HE2);
	}

}


/** Assigens the fixed and moved atoms for ARG
*/
void set_fixed_moved_atoms_side_chains_ARG(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_NE, 
	*atom_CZ, *atom_NH1, *atom_NH2, 
	*atom_HG1, *atom_HG2, *atom_HD1, *atom_HD2, *atom_HE, 
	*atom_1HH1, *atom_2HH1, *atom_1HH2, 
	*atom_2HH2 , *atom_HB1, *atom_HB2;
	int exists_atom_HG1, exists_atom_HG2, exists_atom_HD1, 
	exists_atom_HD2, exists_atom_HE, exists_atom_1HH1, exists_atom_2HH1, 
	exists_atom_1HH2, exists_atom_2HH2, exists_atom_HB1, exists_atom_HB2;
	int num_moved;
	int i_af;
	
	exists_atom_HG1 = 0; 
	exists_atom_HG2 = 0; 
	exists_atom_HD1 = 0;
	exists_atom_HD2 = 0;
	exists_atom_HE = 0;
	exists_atom_1HH1 = 0; 
	exists_atom_2HH1 = 0;
	exists_atom_1HH2 = 0; 
	exists_atom_2HH2 = 0;
	exists_atom_HB1 = 0;
	exists_atom_HB2 = 0;

	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HE = Malloc(char, 3);
		atom_1HH1 = Malloc(char, 5);
		atom_2HH1 = Malloc(char, 5);
		atom_1HH2 = Malloc(char, 5);
		atom_2HH2 = Malloc(char, 5);
		atom_HB1 = Malloc(char, 4);
		atom_HB2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE, "HE");
		strcpy(atom_1HH1, "1HH1");
		strcpy(atom_2HH1, "2HH1");
		strcpy(atom_1HH2, "1HH2");
		strcpy(atom_2HH2, "2HH2");
		strcpy(atom_HB1, "HB1");
		strcpy(atom_HB2, "HB2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CA, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 6;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HB1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HB2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HB2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH2 = 1;
		}		

		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH2, &prot->p_topol->numatom);
		if (exists_atom_HB1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HB2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HB2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);

		}
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);

		}
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE, &prot->p_topol->numatom);

		}
		if (exists_atom_1HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_1HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH2, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH2, &prot->p_topol->numatom);

		}

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);
		free(atom_HG1);
		free(atom_HG2);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HE);
		free(atom_1HH1);
		free(atom_2HH1);
		free(atom_1HH2);
		free(atom_2HH2);		
		free(atom_HB1);
		free(atom_HB2);	

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		atom_HE = Malloc(char, 3);
		atom_1HH1 = Malloc(char, 5);
		atom_2HH1 = Malloc(char, 5);
		atom_1HH2 = Malloc(char, 5);
		atom_2HH2 = Malloc(char, 5);
		atom_HG1 = Malloc(char, 4);
		atom_HG2 = Malloc(char, 4);

		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");
		strcpy(atom_HE, "HE");
		strcpy(atom_1HH1, "1HH1");
		strcpy(atom_2HH1, "2HH1");
		strcpy(atom_1HH2, "1HH2");
		strcpy(atom_2HH2, "2HH2");
		strcpy(atom_HG1, "HG1");
		strcpy(atom_HG2, "HG2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CB, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 5;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HG1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HG2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HG2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH2, &prot->p_topol->numatom);
		if (exists_atom_HG1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HG2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HG2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);

		}
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE, &prot->p_topol->numatom);

		}
		if (exists_atom_1HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_1HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH2, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH2, &prot->p_topol->numatom);

		}


		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);
		free(atom_HD1);
		free(atom_HD2);
		free(atom_HE);
		free(atom_1HH1);
		free(atom_2HH1);
		free(atom_1HH2);
		free(atom_2HH2);		
		free(atom_HG1);
		free(atom_HG2);	
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		atom_HE = Malloc(char, 3);
		atom_1HH1 = Malloc(char, 5);
		atom_2HH1 = Malloc(char, 5);
		atom_1HH2 = Malloc(char, 5);
		atom_2HH2 = Malloc(char, 5);
		atom_HD1 = Malloc(char, 4);
		atom_HD2 = Malloc(char, 4);
		
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");
		strcpy(atom_HE, "HE");
		strcpy(atom_1HH1, "1HH1");
		strcpy(atom_2HH1, "2HH1");
		strcpy(atom_1HH2, "1HH2");
		strcpy(atom_2HH2, "2HH2");
		strcpy(atom_HD1, "HD1");
		strcpy(atom_HD2, "HD2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 4;
		if (atom_name_exists_in_resnum(prot->p_atoms,
				res_num, atom_HD1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HD2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HD2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH2 = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH2, &prot->p_topol->numatom);
		if (exists_atom_HD1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD1, &prot->p_topol->numatom);
		}	
		if (exists_atom_HD2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HD2, &prot->p_topol->numatom);
		}		
		if (exists_atom_HE == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE, &prot->p_topol->numatom);

		}
		if (exists_atom_1HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_1HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH2, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH2, &prot->p_topol->numatom);

		}		

		free(atom_CG);
		free(atom_CD);
		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);
		free(atom_HE);
		free(atom_1HH1);
		free(atom_2HH1);
		free(atom_1HH2);
		free(atom_2HH2);		
		free(atom_HD1);
		free(atom_HD2);	

	}else if (*chi == 4){
				
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_HE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		atom_1HH1 = Malloc(char, 5);
		atom_2HH1 = Malloc(char, 5);
		atom_1HH2 = Malloc(char, 5);
		atom_2HH2 = Malloc(char, 5);

		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_HE, "HE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");
		strcpy(atom_1HH1, "1HH1");
		strcpy(atom_2HH1, "2HH1");
		strcpy(atom_1HH2, "1HH2");
		strcpy(atom_2HH2, "2HH2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE, &prot->p_topol->numatom);
		//Building Moved Atoms
		num_moved = 3;
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH1, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH1 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_1HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_1HH2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_2HH2, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_2HH2 = 1;
		}		
		if (atom_name_exists_in_resnum(prot->p_atoms,
			res_num, atom_HE, &prot->p_topol->numatom) == btrue){
			num_moved = num_moved + 1;
			exists_atom_HE = 1;
		}		
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = num_moved;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH2, &prot->p_topol->numatom);
		if (exists_atom_1HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH1 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH1, &prot->p_topol->numatom);

		}
		if (exists_atom_1HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_1HH2, &prot->p_topol->numatom);

		}
		if (exists_atom_2HH2 == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_2HH2, &prot->p_topol->numatom);

		}
		if (exists_atom_HE == 1){
			i_af++;				
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
					res_num, atom_HE, &prot->p_topol->numatom);

		}
		free(atom_CD);
		free(atom_NE);
		free(atom_HE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);
		free(atom_1HH1);
		free(atom_2HH1);
		free(atom_1HH2);
		free(atom_2HH2);		

	}

}