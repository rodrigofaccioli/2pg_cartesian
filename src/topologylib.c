#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>

#include "defines.h"
#include "messages.h"
#include "futil.h"
#include "topology_types.h"
#include "topologyio.h"
#include "topologylib.h"
#include "topology_atom_parameters.h"
#include "pdbatom.h"
#include "string_owner.h"
#include "osutil.h"
#include "randomlib.h"
#include "vector_math.h"

float _compute_diehdral_angle(const own_vector_t *a1,
		const own_vector_t *a2,const own_vector_t *a3,	const own_vector_t *a4){
	/* Computes the dihedral angle	 */
	own_vector_t *P1, *P2, *M, *r1, *r2, *r3;
	double mod_P1, mod_P2;
	double W;

	P1 = Malloc(own_vector_t,1);
	P2 = Malloc(own_vector_t,1);
	M = Malloc(own_vector_t,1);
	r1 = Malloc(own_vector_t,1);
	r2 = Malloc(own_vector_t,1);
	r3 = Malloc(own_vector_t,1);

	//Computing distances
	sub_vector(r1,a1,a2);
	sub_vector(r2,a2,a3);
	sub_vector(r3,a3,a4);

	cross_product(P1,r1,r2);
	cross_product(P2,r2,r3);
	mod_P1 = mod_vector(P1);
	mod_P2 = mod_vector(P2);

	W = acos(scalar_prod(P1,P2)/(mod_P1*mod_P2));
	//Check if is necessary change the signal of W
	cross_product(M,P1,P2);
	if (scalar_prod(M,r2) < 0){
		W = W *(-1);
	}
	//Deallocating variables
	free(P1);
	free(P2);
	free(M);
	free(r1);
	free(r2);
	free(r3);

	return W;

}


float _compute_phi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global){
	/* Receives the residue, the pdbatoms structure and the phi topology
	 * Returns the value of phi angle
	 * When residue is N-Terminal, there is not a phi angle
	 *
	 * get_pdb_atom_from_resnum_atomid function returns a pdb_atom structure
	 * which is a line of pdbatom section from pdb file. That line contain an
	 * own_vector_t structure that is stored in a1, a2, a3 and a4.
	 */
	if (top_global->top_global_dieh_phi[*r].atom_number1 == 0){//N-Terminal
		return 0;
	}
	const own_vector_t *a1, *a2, *a3, *a4;
	const pdb_atom_t *aux_pdb_atom_1, *aux_pdb_atom_2, *aux_pdb_atom_3,
	*aux_pdb_atom_4 = NULL;
	float phi;
	aux_pdb_atom_1 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number1-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number1-1].atom_id,
			&top_global->numatom);
	a1 = &aux_pdb_atom_1->coord;
	aux_pdb_atom_2 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number2-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number2-1].atom_id,
			&top_global->numatom);
	a2 = &aux_pdb_atom_2->coord;
	aux_pdb_atom_3 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number3-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number3-1].atom_id,
			&top_global->numatom);
	a3 = &aux_pdb_atom_3->coord;
	aux_pdb_atom_4 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number4-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_phi[*r].atom_number4-1].atom_id,
			&top_global->numatom);
	a4 = &aux_pdb_atom_4->coord;
	phi = _compute_diehdral_angle(a1,a2,a3,a4);
	return phi;
}

float _compute_psi(const int *r, const pdb_atom_t *pdb_atoms,
		const top_global_t *top_global){
	/* Receives the residue, the pdbatoms structure and the psi topology
	 * Returns the value of psi angle
	 * When residue is C-Terminal, there is not a psi angle
	 *
	 * get_pdb_atom_from_resnum_atomid function returns a pdb_atom structure
	 * which is a line of pdbatom section from pdb file. That line contain an
	 * own_vector_t structure that is stored in a1, a2, a3 and a4.
	 */
	if (top_global->top_global_dieh_psi[*r].atom_number4 == 0){//C-Terminal
		return 0;
	}
	const own_vector_t *a1, *a2, *a3, *a4;
	const pdb_atom_t *aux_pdb_atom_1, *aux_pdb_atom_2, *aux_pdb_atom_3,
	*aux_pdb_atom_4 = NULL;
	float psi;
	aux_pdb_atom_1 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number1-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number1-1].atom_id,
			&top_global->numatom);
	a1 = &aux_pdb_atom_1->coord;
	aux_pdb_atom_2 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number2-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number2-1].atom_id,
			&top_global->numatom);
	a2 = &aux_pdb_atom_2->coord;
	aux_pdb_atom_3 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number3-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number3-1].atom_id,
			&top_global->numatom);
	a3 = &aux_pdb_atom_3->coord;
	aux_pdb_atom_4 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number4-1].res_number,
			&top_global->top_global_atom[top_global->top_global_dieh_psi[*r].atom_number4-1].atom_id,
			&top_global->numatom);
	a4 = &aux_pdb_atom_4->coord;

	psi = _compute_diehdral_angle(a1,a2,a3,a4);
	return psi;
}

float  _compute_side_chains_angles(const int *r, const int *chi,
		const pdb_atom_t *pdb_atoms, const top_global_t *top_global) {
	/* Receives the residue, the index of chi, the pdbatoms structure and
	 * the side_chains topology
	 * The chi index is based on residue. Your value can be obtained in
	 * topol_residues_ff structure at topology_charmm27_parameters.h
	 *
	 * Returns the value of chi angle
	 *
	 * get_pdb_atom_from_resnum_atomid function returns a pdb_atom structure
	 * which is a line of pdbatom section from pdb file. That line contain an
	 * own_vector_t structure that is stored in a1, a2, a3 and a4.
	 */
	const own_vector_t *a1, *a2, *a3, *a4;
	int index;
	const pdb_atom_t *aux_pdb_atom_1, *aux_pdb_atom_2, *aux_pdb_atom_3,
	*aux_pdb_atom_4 = NULL;
	float chi_angle;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(1,r,
			chi, top_global);
	aux_pdb_atom_1 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a1 = &aux_pdb_atom_1->coord;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(2,r,
			chi, top_global);
	aux_pdb_atom_2 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a2 = &aux_pdb_atom_2->coord;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(3,r,
			chi, top_global);
	aux_pdb_atom_3 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a3 = &aux_pdb_atom_3->coord;
	index = _get_atom_index_from_top_global_dihedral_side_chain_t(4,r,
			chi, top_global);
	aux_pdb_atom_4 = get_pdb_atom_from_resnum_atomid(pdb_atoms,
			&top_global->top_global_atom[index].res_number,
			&top_global->top_global_atom[index].atom_id,
			&top_global->numatom);
	a4 = &aux_pdb_atom_4->coord;
	chi_angle = _compute_diehdral_angle(a1,a2,a3,a4);
	return chi_angle;

}

type_aminos_t _get_amino_id_3(char *c){
/*Receives an amino (char) and returns your id*/
		type_aminos_t amino_id;
		char msg[30];
		if ( (strcmp(c,"GLY") == 0) ){
			amino_id =  aGLY;
		}else if ( (strcmp(c,"ARG") == 0) ){
			amino_id =  aARG;
		}else if ( (strcmp(c,"ALA") == 0) ){
			amino_id =  aALA;
		}else if ( (strcmp(c,"SER") == 0) ){
			amino_id =  aSER;
		}else if ( (strcmp(c,"THR") == 0) ){
			amino_id =  aTHR;
		}else if ( (strcmp(c,"CYS") == 0) ){
			amino_id =  aCYS;
		}else if ( (strcmp(c,"VAL") == 0) ){
			amino_id =  aVAL;
		}else if ( (strcmp(c,"LEU") == 0) ){
			amino_id =  aLEU;
		}else if ( (strcmp(c,"ILE") == 0) ){
			amino_id =  aILE;
		}else if ( (strcmp(c,"MET") == 0) ){
			amino_id =  aMET;
		}else if ( (strcmp(c,"PRO") == 0) ){
			amino_id =  aPRO;
		}else if ( (strcmp(c,"PHE") == 0) ){
			amino_id =  aPHE;
		}else if ( (strcmp(c,"TYR") == 0) ){
			amino_id =  aTYR;
		}else if ( (strcmp(c,"TRP") == 0) ){
			amino_id =  aTRP;
		}else if ( (strcmp(c,"ASP") == 0) ){
			amino_id =  aASP;
		}else if ( (strcmp(c,"GLU") == 0) ){
			amino_id =  aGLU;
		}else if ( (strcmp(c,"ASN") == 0) ){
			amino_id =  aASN;
		}else if ( (strcmp(c,"GLN") == 0) ){
			amino_id =  aGLN;
		}else if ( (strcmp(c,"HIS") == 0) ){
			amino_id =  aHIS;
		}else if ( (strcmp(c,"LYS") == 0) ){
			amino_id =  aLYS;
		}else{
			if (strcmp(c,"") == 0){
				sprintf(msg,"Amino not found, because amino variable is empty. Check it. \n");
			}else{
				sprintf(msg,"%s Amino not found. Check it. \n",c);
			}
			fatal_error(msg);
		}
		return amino_id;
}


type_aminos_t _get_amino_id_1(char c){
/*receives an amino (char) and returns its id*/
	    type_aminos_t amino_id;
		switch (c)   {
			case 'A': //alanina
				amino_id =  aALA;
				break;
			case 'V': //valina
				amino_id =  aVAL;
				break;
			case 'F': //fenilalanina
				amino_id =  aPHE;
				break;
			case 'P': //prolina
				amino_id =  aPRO;
				break;
			case 'L': //leucina
				amino_id =  aLEU;
				break;
			case 'I': //iso leucina
				amino_id =  aILE;
				break;
			case 'R': //argenina
				amino_id =  aARG;
				break;
			case 'D': //acido aspartico
				amino_id =  aASP;
				break;
			case 'E': //acido glutamico
				amino_id =  aGLU;
				break;
			case 'S': //serina
				amino_id =  aSER;
				break;
			case 'T': //treonina
				amino_id =  aTHR;
				break;
			case 'C': //cisteina
				amino_id =  aCYS;
				break;
			case 'N': //asparagina
				amino_id =  aASN;
				break;
			case 'Q': //glutanima
				amino_id =  aGLN;
				break;
			case 'H': //histidina
				amino_id =  aHIS;
				break;
			case 'K': //lisina
				amino_id =  aLYS;
				break;
			case 'Y': //tirosina
				amino_id =  aTYR;
				break;
			case 'M': //metionina
				amino_id =  aMET;
				break;
			case 'W': //triptofano
				amino_id =  aTRP;
				break;
			case 'G': //glicina
				amino_id =  aGLY;
				break;
			default:
				amino_id = aNR;
				char mens [] = "Amino value is not known. Please, check your file.";
				fatal_error(mens);
		}
		return amino_id;
}

void _type_of_diedhral_angle2str(char *str,
		const type_dihedral_angles_t *type_dihedral){
	/* Converts type_dihedral to string
	 * type_dihedral_angles_t is defined at enums.h
	*/
	if ( *type_dihedral == angl_phi){
		strcpy(str,"PHI");
	}else if (*type_dihedral == angl_psi){
		strcpy(str,"PSI");
	}else if (*type_dihedral == angl_chi1){
		strcpy(str,"CHI1");
	}else if (*type_dihedral == angl_chi2){
		strcpy(str,"CHI2");
	}else if (*type_dihedral == angl_chi3){
		strcpy(str,"CHI3");
	}else if (*type_dihedral == angl_chi4){
		strcpy(str,"CHI4");
	}else if (*type_dihedral == angl_typ_dieh_1){
		strcpy(str,"W1");
	}else if (*type_dihedral == angl_typ_dieh_2){
		strcpy(str,"W2");
	}else if (*type_dihedral == angl_typ_dieh_3){
		strcpy(str,"W3");
	}else if (*type_dihedral == angl_typ_dieh_180){
		strcpy(str,"dieh_180");
	}else if (*type_dihedral == angl_typ_dieh_0){
		strcpy(str,"dieh_0");
	}else if (*type_dihedral == angl_typ_dieh_90){
		strcpy(str,"dieh_90");
	}else if (*type_dihedral == angl_typ_dieh_117_){
		strcpy(str,"dieh_117_");
	}else if (*type_dihedral == angl_typ_dieh_trans_120){
		strcpy(str,"dieh_trans_120");
	}else if (*type_dihedral == angl_typ_dieh_trans_240){
		strcpy(str,"dieh_trans_240");
	}else if (*type_dihedral == angl_typ_dieh_trans_123_){
		strcpy(str,"angl_typ_dieh_trans_123_");
	}else{
		fatal_error("Type of dihedral not found. Please look at type_dihedral_angles_t enum.\n");
	}
}

static int _get_atom_index_from_top_global_dihedral_side_chain_t(
		int atm_opt, const int *r, 	const int *chi,
		const top_global_t *top_global){
	/*
	 * Return the index of atom which is used to compute side chains. This
	 * function is used in _compute_side_chains_angles function.
	 * Must be (*r+1) because r is started with 0 and it is employed a lot of
	 * functions. Please look at _pdbatoms2protein function.
	 * Must be (*chi+1) because chi is started with 0. It is an index. Please
	 * look at _pdbatoms2protein function.
	 * The value -1 because index is atom number - 1
	 */
	int aux = -1;
	for (int i =0; i < top_global->number_protein_side_chains; i++){
		if (top_global->top_global_dieh_side_chains[i].res_number == (*r+1)){
			if (top_global->top_global_dieh_side_chains[i].chi == (*chi+1)){
				//Returns based on atom option
				if (atm_opt == 1){
					aux =  top_global->top_global_dieh_side_chains[i].atom_number1 - 1;
				}else if (atm_opt == 2){
					aux = top_global->top_global_dieh_side_chains[i].atom_number2 - 1;
				}else if (atm_opt == 3){
					aux = top_global->top_global_dieh_side_chains[i].atom_number3 - 1;
				}else if (atm_opt == 4){
					aux = top_global->top_global_dieh_side_chains[i].atom_number4 - 1;
				}
			}
		}
	}
	if (aux == -1){
		fatal_error("Index not found in _get_atom_index_from_top_global_dihedral_side_chain_t \n");
	}
	return aux;
}