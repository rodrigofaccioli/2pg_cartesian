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
#include "topology.h"
#include "topologyio.h"
#include "rotation.h"
#include "math_owner.h"
#include "vector_math.h"
#include "solution.h"

float compute_diehdral_angle(const own_vector_t *a1,
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

/* compute phi of residue
*/
float compute_phi_residue(pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top){	
	if (*res_num > 1){
		char *atmC_ = NULL;
		char *atmN = NULL;
		char *atmCA = NULL;
		char *atmC = NULL;
		int res_num_ant;		
		float phi;
		const own_vector_t *a1,*a2,*a3,*a4;

		atmC_ = (char*)malloc(sizeof(char)*2);
		atmN = (char*)malloc(sizeof(char)*2);
		atmCA = (char*)malloc(sizeof(char)*3);
		atmC = (char*)malloc(sizeof(char)*2);
		strcpy(atmC_, "C");
		strcpy(atmN, "N");
		strcpy(atmCA, "CA");
		strcpy(atmC, "C");

		res_num_ant = *res_num - 1;
		a1 = get_pdb_atom_coordinates(prot,&res_num_ant, atmC_, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmN, &top->numatom);
		a3 = get_pdb_atom_coordinates(prot,res_num, atmCA, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmC, &top->numatom);

		phi = compute_diehdral_angle(a1, a2, a3, a4);

		free(atmC_);
		free(atmN);
		free(atmCA);
		free(atmC);

		return phi;
	}
	return 0;
}

/* compute psi of residue
*/
float compute_psi_residue(pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top){	
	if (*res_num < (top->numres-1) ){
		char *atmN = NULL;
		char *atmCA = NULL;
		char *atmC = NULL;
		char *atmN_plus = NULL;

		int res_num_plus;		
		float psi;
		const own_vector_t *a1,*a2,*a3,*a4;
		
		atmN = (char*)malloc(sizeof(char)*2);
		atmCA = (char*)malloc(sizeof(char)*3);
		atmC = (char*)malloc(sizeof(char)*2);
		atmN_plus = (char*)malloc(sizeof(char)*2);		
		strcpy(atmN, "N");
		strcpy(atmCA, "CA");
		strcpy(atmC, "C");
		strcpy(atmN_plus, "N");
		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmN, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmCA, &top->numatom);
		a3 = get_pdb_atom_coordinates(prot,res_num, atmC, &top->numatom);
		res_num_plus = *res_num + 1;
		a4 = get_pdb_atom_coordinates(prot,&res_num_plus, atmN_plus, &top->numatom);
		
		psi = compute_diehdral_angle(a1, a2, a3, a4);
		
		free(atmN);
		free(atmCA);
		free(atmC);
		free(atmN_plus);

		return psi;
	}
	return 0;
}

/* compute omega of residue
*/
float compute_omega_residue(pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top){	
	if (*res_num < (top->numres-1) ){		
		char *atmCA = NULL;
		char *atmC = NULL;
		char *atmCA_plus = NULL;
		char *atmN_plus = NULL;

		int res_num_plus;		
		float omega;
		const own_vector_t *a1,*a2,*a3,*a4;
				
		atmCA = (char*)malloc(sizeof(char)*3);
		atmC = (char*)malloc(sizeof(char)*2);
		atmN_plus = (char*)malloc(sizeof(char)*2);		
		atmCA_plus = (char*)malloc(sizeof(char)*3);

		strcpy(atmCA, "CA");
		strcpy(atmC, "C");		
		strcpy(atmN_plus, "N");
		strcpy(atmCA_plus, "CA");
		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmCA, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmC, &top->numatom);		
		res_num_plus = *res_num + 1;
		a3 = get_pdb_atom_coordinates(prot, &res_num_plus, atmN_plus, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,&res_num_plus, atmCA_plus, &top->numatom);
		
		omega = compute_diehdral_angle(a1, a2, a3, a4);
				
		free(atmCA);
		free(atmC);
		free(atmN_plus);
		free(atmCA_plus);

		return omega;
	}
	return 0;
}

/* compute chi of residue
* chi array that contains each chi angles of residues
* num_chi stores how many chi exists in residue
*/
void compute_chi_residue(float *chi, int *num_chi, pdb_atom_t *prot, 
	const int *res_num, const top_global_t *top){	
	
	char *atmA1 = NULL;
	char *atmA2 = NULL;
	char *atmA3 = NULL;
	char *atmA4 = NULL;
	char *res_name = NULL;
	const own_vector_t *a1,*a2,*a3,*a4;

	*num_chi = -1;
	//Obtaing residue name
	res_name = (char*)malloc(sizeof(char)*4);
	get_res_name_from_res_num(res_name,res_num, prot, &top->numatom);

	atmA1 = (char*)malloc(sizeof(char)*5);
	atmA2 = (char*)malloc(sizeof(char)*5);
	atmA3 = (char*)malloc(sizeof(char)*5);		
	atmA4 = (char*)malloc(sizeof(char)*5);

	if ( strcmp(res_name, "ALA") == 0 ){
		*num_chi = 1;
		strcpy(atmA1, "CB");
		strcpy(atmA2, "C");		
		strcpy(atmA3, "CA");
		strcpy(atmA4, "N");
		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);

		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "ARG") == 0 ){
		*num_chi = 4;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi3		
		strcpy(atmA1, "CB");
		strcpy(atmA2, "CG");		
		strcpy(atmA3, "CD");
		strcpy(atmA4, "NE");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[2] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi4			
		strcpy(atmA1, "CG");
		strcpy(atmA2, "CD");		
		strcpy(atmA3, "NE");
		strcpy(atmA4, "CZ");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[3] = compute_diehdral_angle(a1, a2, a3, a4);

	}else if ( strcmp(res_name, "LYS") == 0 ){
		*num_chi = 4;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi3		
		strcpy(atmA1, "CB");
		strcpy(atmA2, "CG");		
		strcpy(atmA3, "CD");
		strcpy(atmA4, "CE");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[2] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi4		
		strcpy(atmA1, "CG");
		strcpy(atmA2, "CD");		
		strcpy(atmA3, "CE");
		strcpy(atmA4, "NZ");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[3] = compute_diehdral_angle(a1, a2, a3, a4);

	}else if ( strcmp(res_name, "ASP") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "OD1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "GLU") == 0 ){
		*num_chi = 3;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi3				
		strcpy(atmA1, "CB");
		strcpy(atmA2, "CG");		
		strcpy(atmA3, "CD");
		strcpy(atmA4, "OE1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[2] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "GLN") == 0 ){
		*num_chi = 3;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi3		
		strcpy(atmA1, "CB");
		strcpy(atmA2, "CG");		
		strcpy(atmA3, "CD");
		strcpy(atmA4, "OE1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[2] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "ASN") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "OD1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "HIS") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "ND1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "SER") == 0 ){
		*num_chi = 1;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "OG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "THR") == 0 ){
		*num_chi = 1;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "OG1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "TYR") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "CYS") == 0 ){
		*num_chi = 1;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "SG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "MET") == 0 ){
		*num_chi = 3;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "SD");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi3		
		strcpy(atmA1, "CB");
		strcpy(atmA2, "CG");		
		strcpy(atmA3, "SD");
		strcpy(atmA4, "CE");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[2] = compute_diehdral_angle(a1, a2, a3, a4);

	}else if ( strcmp(res_name, "TRP") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "ILE") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG1");
		strcpy(atmA4, "CD");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "LEU") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "PHE") == 0 ){
		*num_chi = 2;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);

		//chi2		
		strcpy(atmA1, "CA");
		strcpy(atmA2, "CB");		
		strcpy(atmA3, "CG");
		strcpy(atmA4, "CD1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[1] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "VAL") == 0 ){
		*num_chi = 1;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG1");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "PRO") == 0 ){
		*num_chi = 1;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "ALA") == 0 ){
		*num_chi = 1;
		//chi1		
		strcpy(atmA1, "N");
		strcpy(atmA2, "CA");		
		strcpy(atmA3, "CB");
		strcpy(atmA4, "CG");		
		a1 = get_pdb_atom_coordinates(prot,res_num, atmA1, &top->numatom);
		a2 = get_pdb_atom_coordinates(prot,res_num, atmA2, &top->numatom);			
		a3 = get_pdb_atom_coordinates(prot,res_num, atmA3, &top->numatom);
		a4 = get_pdb_atom_coordinates(prot,res_num, atmA4, &top->numatom);
		chi[0] = compute_diehdral_angle(a1, a2, a3, a4);
	}else if ( strcmp(res_name, "GLY") == 0 ){
		*num_chi = 0;
	}

	free(atmA1);
	free(atmA2);
	free(atmA3);
	free(atmA4);
	free(res_name);
		
	
}

int main(int argc, char *argv[]){
	input_parameters_t *in_para;
	in_para = (input_parameters_t *)malloc(sizeof(input_parameters_t));
	display_msg("Reading the configure file \n");
	load_parameters_from_file(in_para,argv[1]);

    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    
    protein_t *population_p; // main population

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

	//Allocating PDB ATOMS of population_p
    population_p = allocateProtein(&in_para->size_population);    

    //Loading initial population and allocating atom and topology
    load_initial_population_file(population_p, &in_para->size_population, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);

    int res_num, n_chi, ind;
    float phi;
    float psi;
    float omega;
    float *chi = (float*)malloc(sizeof(float)*MAX_CHI);
    int chi_number;

    for (ind = 0; ind < in_para->size_population; ind++){
    	printf("Individuo %i\n", ind);    	
    	for (res_num = 1; res_num <= population_p[ind].p_topol->numres; res_num++){
    		printf("phi\tpsi\tomega do residuo %i \n", res_num);
    		phi = compute_phi_residue(population_p[ind].p_atoms, &res_num, population_p[ind].p_topol);
    		psi = compute_psi_residue(population_p[ind].p_atoms, &res_num, population_p[ind].p_topol);
	    	omega = compute_omega_residue(population_p[ind].p_atoms, &res_num, population_p[ind].p_topol);
    		printf("%f\t%f\t%f\n", phi, psi, omega);    	
    		compute_chi_residue(chi, &chi_number, population_p[ind].p_atoms, &res_num, population_p[ind].p_topol);
			printf("Chi do residuo %i\n", res_num);	
			for (n_chi = 0; n_chi < chi_number; n_chi++){
				printf("%f\n", chi[n_chi]);
			}    	
	    }    
    }
    free(chi);
	deAllocateload_parameters(in_para);

	display_msg("Done 2PG Crossover \n");

	return 0;
}