#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "topologylib.h"
#include "messages.h"
#include "defines.h"
#include "pdbatom.h"

int get_atom_index_by_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname, const int *num_atom){
	const pdb_atom_t * aux;

	aux = search_pdb_atom_from_resnum_atom_name(atoms, res_num, atomname, num_atom);
	if (aux == NULL){
		fatal_error("Atom not found at get_atom_index_by_resnum_atom_name\n");
	}
	return aux->atmnumber;
}


/** Assigens the fixed and moved atoms for SER
*/
void set_fixed_moved_atoms_side_chains_SER(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_OG;
	int i_af;
	
	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_OG = Malloc(char, 3);

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_OG, "OG");

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
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_OG, &prot->p_topol->numatom);

	free(atom_CA);
	free(atom_CB);
	free(atom_OG);
}

/** Assigens the fixed and moved atoms for CYS
*/
void set_fixed_moved_atoms_side_chains_CYS(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_SG;
	int i_af;
	
	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_SG = Malloc(char, 3);

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_SG, "SG");

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
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_SG, &prot->p_topol->numatom);

	free(atom_CA);
	free(atom_CB);
	free(atom_SG);
}

/** Assigens the fixed and moved atoms for THR
*/
void set_fixed_moved_atoms_side_chains_THR(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_OG1, *atom_CG2;
	int i_af;
	
	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_OG1 = Malloc(char, 4);
	atom_CG2 = Malloc(char, 4);

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_OG1, "OG1");
	strcpy(atom_CG2, "CG2");

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
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_OG1, &prot->p_topol->numatom);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CG2, &prot->p_topol->numatom);


	free(atom_CA);
	free(atom_CB);
	free(atom_OG1);
	free(atom_CG2);
}

/** Assigens the fixed and moved atoms for VAL
*/
void set_fixed_moved_atoms_side_chains_VAL(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG1, *atom_CG2;
	int i_af;
	
	atom_CA = Malloc(char, 3);
	atom_CB = Malloc(char, 3);
	atom_CG1 = Malloc(char, 4);
	atom_CG2 = Malloc(char, 4);

	strcpy(atom_CA, "CA");
	strcpy(atom_CB, "CB");
	strcpy(atom_CG1, "CG1");
	strcpy(atom_CG2, "CG2");

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
	i_af = -1;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CG1, &prot->p_topol->numatom);
	i_af++;				
	prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
			res_num, atom_CG2, &prot->p_topol->numatom);


	free(atom_CA);
	free(atom_CB);
	free(atom_CG1);
	free(atom_CG2);
}

/** Assigens the fixed and moved atoms for ASP
*/
void set_fixed_moved_atoms_side_chains_ASP(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_OD1, *atom_OD2;
	int i_af;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_OD2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_OD2, "OD2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_OD2);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_OD2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_OD2, "OD2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD2, &prot->p_topol->numatom);
		
		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_OD2);

	}

}

/** Assigens the fixed and moved atoms for ASN
*/
void set_fixed_moved_atoms_side_chains_ASN(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_OD1, *atom_ND2;
	int i_af;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_ND2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_ND2, "ND2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_ND2);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_OD1 = Malloc(char, 4);
		atom_ND2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_OD1, "OD1");
		strcpy(atom_ND2, "ND2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_ND2, &prot->p_topol->numatom);
		
		free(atom_CB);
		free(atom_CG);
		free(atom_OD1);
		free(atom_ND2);

	}

}

/** Assigens the fixed and moved atoms for LEU
*/
void set_fixed_moved_atoms_side_chains_LEU(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2;
	int i_af;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1, "CD1");
		strcpy(atom_CD2, "CD2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1, "CD1");
		strcpy(atom_CD2, "CD2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD2, &prot->p_topol->numatom);
		
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);

	}

}


/** Assigens the fixed and moved atoms for PRO
*/
void set_fixed_moved_atoms_side_chains_PRO(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD;
	int i_af;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);		

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");	

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CG, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);		

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");		

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CD, &prot->p_topol->numatom);
		
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);		
	}

}

/** Assigens the fixed and moved atoms for PHE
*/
void set_fixed_moved_atoms_side_chains_PHE(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2, 
	*atom_CE1, *atom_CE2, *atom_CZ;
	int i_af;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ = Malloc(char, 3);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ,"CZ" );

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 6;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ = Malloc(char, 3);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ,"CZ" );

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 5;
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

		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);

		
	}

}


/** Assigens the fixed and moved atoms for HIS
*/
void set_fixed_moved_atoms_side_chains_HIS(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_ND1, *atom_CD2, *atom_CE1, *atom_NE2;
	int i_af;
	
	if (*chi == 1){
		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_ND1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_ND1, "ND1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_CE1, "CE1");
		strcpy(atom_NE2, "NE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 5;
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


		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_ND1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_NE2);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_ND1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_ND1, "ND1");
		strcpy(atom_CD2, "CD2");
		strcpy(atom_CE1, "CE1");
		strcpy(atom_NE2, "NE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 5;
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

		free(atom_CB);
		free(atom_CG);
		free(atom_ND1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_NE2);
	}

}


/** Assigens the fixed and moved atoms for TYR
*/
void set_fixed_moved_atoms_side_chains_TYR(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2, 
	*atom_CE1, *atom_CE2, *atom_CZ, *atom_OH;
	int i_af;
	
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

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ,"CZ" );
		strcpy(atom_OH,"OH" );

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 7;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);
		free(atom_OH);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD1 = Malloc(char, 4);
		atom_CD2 = Malloc(char, 4);
		atom_CE1 = Malloc(char, 4);
		atom_CE2 = Malloc(char, 4);
		atom_CZ = Malloc(char, 3);
		atom_OH = Malloc(char, 3);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD1,"CD1" );
		strcpy(atom_CD2,"CD2" );
		strcpy(atom_CE1,"CE1" );
		strcpy(atom_CE2,"CE2" );
		strcpy(atom_CZ,"CZ" );
		strcpy(atom_OH,"OH" );

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 6;
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

		free(atom_CB);
		free(atom_CG);
		free(atom_CD1);
		free(atom_CD2);
		free(atom_CE1);
		free(atom_CE2);
		free(atom_CZ);
		free(atom_OH);
		
	}

}


/** Assigens the fixed and moved atoms for TRP
*/
void set_fixed_moved_atoms_side_chains_TRP(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD1, *atom_CD2, *atom_NE1, *atom_CE3, 
	*atom_CE2, *atom_CZ3, *atom_CZ2, *atom_CH2;
	int i_af;
	
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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 9;
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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 8;
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
	}

}


/** Assigens the fixed and moved atoms for MET
*/
void set_fixed_moved_atoms_side_chains_MET(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_SD, *atom_CE;
	int i_af;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_SD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_SD, "SD");
		strcpy(atom_CE, "CE");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_SD);
		free(atom_CE);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_SD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);

		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_SD, "SD");
		strcpy(atom_CE, "CE");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_SD, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
	
		free(atom_CB);
		free(atom_CG);
		free(atom_SD);
		free(atom_CE);
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_SD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		
		strcpy(atom_CG, "CG");
		strcpy(atom_SD, "SD");
		strcpy(atom_CE, "CE");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
			
		free(atom_CG);
		free(atom_SD);
		free(atom_CE);	
	}

}

/** Assigens the fixed and moved atoms for GLN
*/
void set_fixed_moved_atoms_side_chains_GLN(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_OE1, *atom_NE2;
	int i_af;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_NE2, "NE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 4;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_NE2);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_NE2, "NE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_NE2);
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_NE2 = Malloc(char, 4);
				
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_NE2, "NE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE2, &prot->p_topol->numatom);

		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_NE2);

	}

}


/** Assigens the fixed and moved atoms for GLU
*/
void set_fixed_moved_atoms_side_chains_GLU(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_OE1, *atom_OE2;
	int i_af;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_OE2 = Malloc(char, 4);

		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_OE2, "OE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 4;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_OE2);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_OE2 = Malloc(char, 4);
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_OE2, "OE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_OE2);
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_OE1 = Malloc(char, 4);
		atom_OE2 = Malloc(char, 4);
				
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_OE1, "OE1");
		strcpy(atom_OE2, "OE2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_OE2, &prot->p_topol->numatom);

		free(atom_CG);
		free(atom_CD);
		free(atom_OE1);
		free(atom_OE2);

	}

}


/** Assigens the fixed and moved atoms for LYS
*/
void set_fixed_moved_atoms_side_chains_LYS(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_CE, *atom_NZ;
	int i_af;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char, 3);
		
		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 4;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char, 3);
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char,3);
				
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CE, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NZ, &prot->p_topol->numatom);

		free(atom_CG);
		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);

	}else if (*chi == 4){
				
		atom_CD = Malloc(char, 3);
		atom_CE = Malloc(char, 3);
		atom_NZ = Malloc(char,3);
						
		strcpy(atom_CD, "CD");
		strcpy(atom_CE, "CE");
		strcpy(atom_NZ, "NZ");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NZ, &prot->p_topol->numatom);

		free(atom_CD);
		free(atom_CE);
		free(atom_NZ);

	}

}


/** Assigens the fixed and moved atoms for ARG
*/
void set_fixed_moved_atoms_side_chains_ARG(protein_t *prot, 
	const int *res_num, const int *chi){

	char *atom_CA, *atom_CB, *atom_CG, *atom_CD, *atom_NE, 
	*atom_CZ, *atom_NH1, *atom_NH2;
	int i_af;
	
	if (*chi == 1){

		atom_CA = Malloc(char, 3);
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		
		strcpy(atom_CA, "CA");
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 6;
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

		free(atom_CA);
		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);

	}else if (*chi == 2){
		
		atom_CB = Malloc(char, 3);
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		
		strcpy(atom_CB, "CB");
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 5;
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

		free(atom_CB);
		free(atom_CG);
		free(atom_CD);
		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);
		
	}else if (*chi == 3){
				
		atom_CG = Malloc(char, 3);
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		
		strcpy(atom_CG, "CG");
		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 4;
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

		free(atom_CG);
		free(atom_CD);
		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);

	}else if (*chi == 4){
				
		atom_CD = Malloc(char, 3);
		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		
		strcpy(atom_CD, "CD");
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");

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
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 3;
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

		free(atom_CD);
		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);

	}else if (*chi == 5){

		atom_NE = Malloc(char, 3);
		atom_CZ = Malloc(char, 3);
		atom_NH1 = Malloc(char, 4);
		atom_NH2 = Malloc(char, 4);
		
		strcpy(atom_NE, "NE");
		strcpy(atom_CZ, "CZ");
		strcpy(atom_NH1, "NH1");
		strcpy(atom_NH2, "NH2");

		//Building Fixed Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_fixed);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NE, &prot->p_topol->numatom);
		i_af++;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].fixed_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_CZ, &prot->p_topol->numatom);
		//Building Moved Atoms
		i_af = -1;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved = 2;
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms = Malloc(int, 
			prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].num_moved);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH1, &prot->p_topol->numatom);
		i_af++;				
		prot->p_topol->side_chains[*res_num-1].atoms_chi[*chi-1].moved_atoms[i_af] = get_atom_index_by_resnum_atom_name(prot->p_atoms,
				res_num, atom_NH2, &prot->p_topol->numatom);

		free(atom_NE);
		free(atom_CZ);
		free(atom_NH1);
		free(atom_NH2);

	}

}