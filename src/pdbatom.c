#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "pdbatom.h"
#include "defines.h"
#include "string_owner.h"
#include "messages.h"

_2PG_CARTESIAN_EXPORT
pdb_atom_t** allocate_Population_pdb(const int *inPopSize, const int *numatom_by_model){
	/* Allocate a population of pdb_atom_t
	*/
	 pdb_atom_t **population;
     population = Malloc(pdb_atom_t*,*inPopSize);
     for(int i=0;i<*inPopSize;i++){
    	population[i] = allocate_pdbatom(&numatom_by_model[i]);
 	}
    return population;
}

_2PG_CARTESIAN_EXPORT
void desAllocate_Population_pdb(pdb_atom_t** pdbatoms, const int *inPopSize){
	for (int i = 0; i < *inPopSize; i++){
		desAllocate_pdbatom(pdbatoms[i]);
	}	
}

_2PG_CARTESIAN_EXPORT
pdb_atom_t * allocate_pdbatom(const int *numatom){
	pdb_atom_t *pdbatom_aux = NULL;
	pdbatom_aux = Malloc(pdb_atom_t,*numatom);
	return pdbatom_aux;
}

_2PG_CARTESIAN_EXPORT
void desAllocate_pdbatom(pdb_atom_t *pdbatoms){
	free(pdbatoms);
	pdbatoms = NULL;
}

_2PG_CARTESIAN_EXPORT
void set_pdb_atom_coordinates(pdb_atom_t *pdbatom,	char *atmname, char *resname,
		const char *chain_name, const int *resnum, 	const float *x,
		const float *y, const float *z, const int *index){
	set_pdb_atom_generic_information(pdbatom, atmname, resname, chain_name,
			resnum, index);
	pdbatom->coord.x = *x;
	pdbatom->coord.y = *y;
	pdbatom->coord.z = *z;
}

_2PG_CARTESIAN_EXPORT
void set_pdb_atom_generic_information(pdb_atom_t *pdbatom,
		char *atmname, char *resname,
		const char *chain_name, const int *resnum,  const int *index){
	/* Receives values without coordinates */
	trim(atmname);
	trim(resname);
	strcpy(pdbatom->atmname, atmname);
	strcpy(pdbatom->resname, resname);
	pdbatom->resnum  = *resnum;
	pdbatom->atomid = get_atomid_from_atom_name(atmname);
	pdbatom->atmnumber = *index;
}

void get_atom_name_from_atomid(char *atomname, const type_atoms_t *atomid){
	if (*atomid == atmC){
		strcpy(atomname,"C");
	}else if (*atomid == atmCA){
		strcpy(atomname,"CA");
	}else if (*atomid == atmN){
		strcpy(atomname,"N");
	}else if (*atomid == atmCB){
		strcpy(atomname,"CB");
	}else if (*atomid == atmCG){
		strcpy(atomname,"CG");
	}else if (*atomid == atmCD){
		strcpy(atomname,"CD");
	}else if (*atomid == atmNE){
		strcpy(atomname,"NE");
	}else if (*atomid == atmCZ){
		strcpy(atomname,"CZ");
	}else if (*atomid == atmOD1){
		strcpy(atomname,"OD1");
	}else if (*atomid == atmSG){
		strcpy(atomname,"SG");
	}else if (*atomid == atmOE1){
		strcpy(atomname,"OE1");
	}else if (*atomid == atmCG1){
		strcpy(atomname,"CG1");
	}else if (*atomid == atmCD1){
		strcpy(atomname,"CD1");
	}else if (*atomid == atmCE){
		strcpy(atomname,"CE");
	}else if (*atomid == atmNZ){
		strcpy(atomname,"NZ");
	}else if (*atomid == atmOG){
		strcpy(atomname,"OG");
	}else if (*atomid == atmSD){
		strcpy(atomname,"SD");
	}else if (*atomid == atmOG1){
		strcpy(atomname,"OG1");
	}else if (*atomid == atmND1){
		strcpy(atomname,"ND1");
	}else{
		char msg[300];
		sprintf(msg,"Not found out atom, check get_atom_name_from_atomid function.\n");
		fatal_error(msg);
	}

}

type_atoms_t get_atomid_from_atom_name(const  char *__atmname){
	/*Receives atom name
	 * Returns its atom id
	 */
	char msg[300];
	if (is_equal(__atmname,"C")){
		return atmC;
	}else if (is_equal(__atmname,"CA")){
		return atmCA;
	}else if (is_equal(__atmname,"CB")){
		return atmCB;
	}else if (is_equal(__atmname,"CD")){
		return atmCD;
	}else if (is_equal(__atmname,"CD1")){
		return atmCD1;
	}else if (is_equal(__atmname,"CD2")){
		return atmCD2;
	}else if (is_equal(__atmname,"CE")){
		return atmCE;
	}else if (is_equal(__atmname,"CE1")){
		return atmCE1;
	}else if (is_equal(__atmname,"CE2")){
		return atmCE2;
	}else if (is_equal(__atmname,"CE3")){
		return atmCE3;
	}else if (is_equal(__atmname,"CG")){
		return atmCG;
	}else if (is_equal(__atmname,"CG1")){
		return atmCG1;
	}else if (is_equal(__atmname,"CG2")){
		return atmCG2;
	}else if (is_equal(__atmname,"CH2")){
		return atmCH2;
	}else if (is_equal(__atmname,"CZ")){
		return atmCZ;
	}else if (is_equal(__atmname,"N")){
		return atmN;
	}else if (is_equal(__atmname,"ND1")){
		return atmND1;
	}else if (is_equal(__atmname,"ND2")){
		return atmND2;
	}else if (is_equal(__atmname,"NE")){
		return atmNE;
	}else if (is_equal(__atmname,"NE1")){
		return atmNE1;
	}else if (is_equal(__atmname,"NE2")){
		return atmNE2;
	}else if (is_equal(__atmname,"NH1")){
		return atmNH1;
	}else if (is_equal(__atmname,"NH2")){
		return atmNH2;
	}else if (is_equal(__atmname,"NZ")){
		return atmNZ;
	}else if (is_equal(__atmname,"CZ")){
		return atmCZ;
	}else if (is_equal(__atmname,"CZ2")){
		return atmCZ2;
	}else if (is_equal(__atmname,"CZ3")){
		return atmCZ3;
	}else if (is_equal(__atmname,"O")){
		return atmO;
	}else if (is_equal(__atmname,"OD1")){
		return atmOD1;
	}else if (is_equal(__atmname,"OD2")){
		return atmOD2;
	}else if (is_equal(__atmname,"H")){
		return atmH;
	}else if (is_equal(__atmname,"HA")){
		return atmHA;
	}else if (is_equal(__atmname,"HA1")){
		return atmHA1;
	}else if (is_equal(__atmname,"HA2")){
		return atmHA2;
	}else if (is_equal(__atmname,"HA3")){
		return atmHA3;
	}else if (is_equal(__atmname,"HE")){
		return atmHE;
	}else if (is_equal(__atmname,"HE1")){
		return atmHE1;
	}else if (is_equal(__atmname,"HE2")){
		return atmHE2;
	}else if (is_equal(__atmname,"HE3")){
		return atmHE3;
	}else if (is_equal(__atmname,"HG")){
		return atmHG;
	}else if (is_equal(__atmname,"HG1")){
		return atmHG1;
	}else if (is_equal(__atmname,"HG2")){
		return atmHG2;
	}else if (is_equal(__atmname,"HG3")){
		return atmHG3;
	}else if (is_equal(__atmname,"HN")){
		return atmHN;
	}else if (is_equal(__atmname,"H1")){
		return atmH1;
	}else if (is_equal(__atmname,"H2")){
		return atmH2;
	}else if (is_equal(__atmname,"H3")){
		return atmH3;
	}else if (is_equal(__atmname,"HB")){
		return atmHB;
	}else if (is_equal(__atmname,"HB1")){
		return atmHB1;
	}else if (is_equal(__atmname,"HB2")){
		return atmHB2;
	}else if (is_equal(__atmname,"HB3")){
		return atmHB3;
	}else if (is_equal(__atmname,"HD1")){
		return atmHD1;
	}else if (is_equal(__atmname,"HD11") ||
			   is_equal(__atmname,"1HD1")){
		return atmHD11;
	}else if (is_equal(__atmname,"HD12") ||
			   is_equal(__atmname,"2HD1")){
		return atmHD12;
	}else if (is_equal(__atmname,"HD13") ||
			   is_equal(__atmname,"3HD1")){
		return atmHD13;
	}else if (is_equal(__atmname,"HD2")){
		return atmHD2;
	}else if (is_equal(__atmname,"HD21") ||
			is_equal(__atmname,"1HD2") ){
		return atmHD21;
	}else if (is_equal(__atmname,"HD22")||
			is_equal(__atmname,"2HD2") 	){
		return atmHD22;
	}else if (is_equal(__atmname,"HD23")||
			is_equal(__atmname,"3HD2") ){
		return atmHD23;
	}else if (is_equal(__atmname,"HD3")){
		return atmHD3;
	}else if ( is_equal(__atmname,"HH11") ||
			   is_equal(__atmname,"1HH1")){
		return atmHH11;
	}else if (is_equal(__atmname,"HH12") ||
			  is_equal(__atmname,"2HH1")){
		return atmHH12;
	}else if (is_equal(__atmname,"HH2")){
		return atmHH2;
	}else if (is_equal(__atmname,"HH21") ||
			  is_equal(__atmname,"1HH2")){
		return atmHH21;
	}else if (is_equal(__atmname,"HH22") ||
			  is_equal(__atmname,"2HH2")){
		return atmHH22;
	}else if (is_equal(__atmname,"HH")){
		return atmHH;
	}else if (is_equal(__atmname,"HZ")){
		return atmHZ;
	}else if (is_equal(__atmname,"HZ1")){
		return atmHZ1;
	}else if (is_equal(__atmname,"HZ2")){
		return atmHZ2;
	}else if (is_equal(__atmname,"HZ3")){
		return atmHZ1;
	}else if (is_equal(__atmname,"OH")){
		return atmOH;
	}else if (is_equal(__atmname,"OT1")){
		return atmOT1;
	}else if (is_equal(__atmname,"OT2")){
		return atmOT2;
	}else if (is_equal(__atmname,"OC1")){
		return atmOC1;
	}else if (is_equal(__atmname,"OC2")){
		return atmOC2;
	}else if (is_equal(__atmname,"OXT")){
		return atmOXT;
	}else if (is_equal(__atmname,"HD21") ||
			  is_equal(__atmname,"1HD2")){
		return atmHD21;
	}else if (is_equal(__atmname,"HD22") ||
			is_equal(__atmname,"2HD2")){
		return atmHD22;
	}else if (is_equal(__atmname,"SG")){
		return atmSG;
	}else if (is_equal(__atmname,"OE1")){
		return atmOE1;
	}else if (is_equal(__atmname,"SD")){
		return atmSD;
	}else if (is_equal(__atmname,"OE2")){
		return atmOE2;
	}else if (is_equal(__atmname,"OG")){
		return atmOG;
	}else if (is_equal(__atmname,"OG1")){
		return atmOG1;
	}else if (is_equal(__atmname,"NE2")){
		return atmNE2;
	}else if (is_equal(__atmname,"HE21")||
			  is_equal(__atmname,"1HE2")){
		return atmHE21;
	}else if (is_equal(__atmname,"HE22")||
			  is_equal(__atmname,"2HE2")){
		return atmHE22;
	}else if (is_equal(__atmname,"HG11")||
			  is_equal(__atmname,"1HG1")){
		return atmHG11;
	}else if (is_equal(__atmname,"HG12")||
			  is_equal(__atmname,"1HG2")){
		return atmHG12;
	}else if (is_equal(__atmname,"HG13")||
			  is_equal(__atmname,"3HG1")){
		return atmHG13;
	}else if (is_equal(__atmname,"HG21")||
			  is_equal(__atmname,"2HG1")){
		return atmHG21;
	}else if (is_equal(__atmname,"HG22")||
			  is_equal(__atmname,"2HG2")){
		return atmHG22;
	}else if (is_equal(__atmname,"HG23")||
			  is_equal(__atmname,"3HG2")){
		return atmHG23;
	}else if (is_equal(__atmname,"HH31")||
			  is_equal(__atmname,"1HH3")){
		return atmHH31;
	}else if (is_equal(__atmname,"HH32")||
			  is_equal(__atmname,"2HH3")){
		return atmHH32;
	}else if (is_equal(__atmname,"HH33")||
			  is_equal(__atmname,"3HH3")){
		return atmHH33;
	}else if (is_equal(__atmname,"CH3")){
		return atmCH3;
	}else{
		sprintf(msg,"Not found out atom %s, check get_atomid_from_atom_name function.\n",__atmname);
		fatal_error(msg);
	}
		
}

static pdb_atom_t * search_pdb_atom_from_resnum_atomid_alow_change(pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom){
	/* Returns pdb_atom_t  based on residue number and atom id. Otherwise,
	 * returns NULL if not found.
	 * It allows to change information about the atom
	 */
	int i;
	for (i = 0; i < *num_atom;i++){
		if (atoms[i].resnum ==  *res_num){
			if (atoms[i].atomid == *atomid){
				return &atoms[i];
			}
		}
	}
	return NULL;
}

const pdb_atom_t * search_pdb_atom_from_resnum_atomid(const pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom){
	/*Returns pdb_atom_t  based on residue number and atom id. Otherwise,
	 * returns NULL if not found.
	 */
	int i;
	for (i = 0; i < *num_atom;i++){
		if (atoms[i].resnum ==  *res_num){
			if (atoms[i].atomid == *atomid){
				return &atoms[i];
			}
		}
	}
	return NULL;
}

const pdb_atom_t * search_pdb_atom_from_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname,	const int *num_atom){
	/*Returns pdb_atom_t  based on residue number and atom name. Otherwise,
	 * returns NULL if not found.
	 */
	int i;
	for (i = 0; i < *num_atom;i++){
		if (atoms[i].resnum ==  *res_num){
			if (strcmp(atoms[i].atmname,atomname) == 0){
				return &atoms[i];
			}
		}
	}
	return NULL;
}

const pdb_atom_t * get_pdb_atom_from_resnum_atomid(const pdb_atom_t *atoms,
		const int *res_num, const type_atoms_t *atomid,	const int *num_atom){
	/* Returns pdb_atom_t  based on residue number and atom id */
	const pdb_atom_t *aux_ret = NULL;
	aux_ret = search_pdb_atom_from_resnum_atomid(atoms, res_num, atomid,
			num_atom);
	if (aux_ret == NULL){
		char msg[300];
		sprintf(msg,"Atom %d not found in residue number %d when was executing get_pdb_atom_from_resnum_atomid function \n",*atomid,*res_num);
		fatal_error(msg);
	}
	return aux_ret;
}

const pdb_atom_t * get_pdb_atom_from_resnum_atom_name(const pdb_atom_t *atoms,
		const int *res_num, const char *atomname,	const int *num_atom){
	/* Returns pdb_atom_t  based on residue number and atom name */
	const pdb_atom_t *aux_ret = NULL;
	aux_ret = search_pdb_atom_from_resnum_atom_name(atoms, res_num, atomname,
			num_atom);
	if (aux_ret == NULL){
		char msg[300];
		sprintf(msg,"Atom %s not found in residue number %d when was executing get_pdb_atom_from_resnum_atom_name function \n",atomname,*res_num);
		fatal_error(msg);
	}
	return aux_ret;
}

/** Based on the number of atom, returns resnum from its last atom which
* represents the number of residues
*/
int get_last_number_res_from_atom(const pdb_atom_t *atoms, const int *num_atom){
	return atoms[*num_atom - 1].resnum;
}

/** Based on the number of atom, returns number of residues 
*/
int get_number_residues_from_atom(const pdb_atom_t *atoms, const int *num_atom){
	int num_res = 0;
	int ref_res;
	int r = 0;	
	while (r < *num_atom){
		ref_res = atoms[r].resnum;
		num_res = num_res + 1;
		while ( (ref_res == atoms[r].resnum) && (r < *num_atom) ){
			r = r + 1;
		}
	}
	return num_res;
}


/** Informs the residue name from residue number
* res_name assined the name of residue
* num_res the number of residue which wants to know the name.
* atoms means the atoms of protein.
* num_atom number of atoms that composes atoms.
*/
void get_res_name_from_res_num(char *res_name, const int *num_res,
		const pdb_atom_t *atoms, const int *num_atom){	
	int a;
	for (a = 0; a < *num_atom; a++){
		if (atoms[a].resnum == *num_res){
			strcpy(res_name,atoms[a].resname);
			break;
		}
	}
}

/** This function changes the resnum when it is incorrect. resnum must be
 * started with 1 for the first residue. See 1VII pdbid. This pdb the first
 * residue is started with 41
 */
void renumerate_residue_number(pdb_atom_t *atoms, const int *num_atom){
	int res_ref, a, res_new_number;
	if (is_residue_number_ok(atoms) == bfalse){
		a = 0;
		res_new_number = 0;
		while (a < *num_atom){
			res_ref = atoms[a].resnum;
			res_new_number = res_new_number + 1;
			while (res_ref == atoms[a].resnum){
				atoms[a].resnum = res_new_number;
				a = a + 1;
			}
		}
	}
}

static boolean_t is_residue_number_ok(pdb_atom_t *atoms){
	/*Checks if residue numbers are fine. Its considered fine
	 * when the resnum from the first atom is 1. Otherwise, is considered
	 * incorrect. Therefore, this functions returns false
	 */
	if (atoms[0].resnum == 1){
		return btrue;
	}else{
		return bfalse;
	}
}

void rename_atom(pdb_atom_t *atoms, const char *name, const char *name_new,
		const int *res_num, const type_atoms_t *atomid,
		const type_atoms_t *atomid_new, const int *num_atom){
	pdb_atom_t * aux;
	aux= search_pdb_atom_from_resnum_atomid_alow_change(atoms, res_num, atomid,
			num_atom);
	if (aux == NULL){
		fatal_error("Atom not found at rename_atom function \n");
	}
	if (strcmp(aux->atmname, name) != 0){
		fatal_error("The atom names do not match at rename_atom function. Check it \n");
	}
	strcpy(aux->atmname,name_new);
	aux->atomid = *atomid_new;
}


_2PG_CARTESIAN_EXPORT
void copy_pdb_atom(pdb_atom_t *dest, const pdb_atom_t *source){
	strcpy(dest->atmname, source->atmname);
	strcpy(dest->resname, source->resname);
	dest->resnum    = source->resnum;
	dest->atomid    = source->atomid;
	dest->atmnumber = source->atmnumber;
	dest->coord.x   = source->coord.x; 
	dest->coord.y   = source->coord.y; 
	dest->coord.z   = source->coord.z; 
}

