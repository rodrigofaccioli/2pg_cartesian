#include <stdlib.h>
#include <math.h>


#include "rotation.h"
#include "vector_types.h"


/** Makes the rotation
* p_atoms represents all atoms
* atmB is fixed atom that was obtained from topology (fixed[0])
* atmC is fixed atom that was obtained from topology (fixed[1])
* atmD is the atom which will be rotated
* dih is the value of rotated angle 
*/
static void rotation_by_angle_dih(pdb_atom_t *p_atoms, const int *atmB, const int *atmC, 
	const int *atmD, const float *dih){
	own_vector_t orign, bc, u;

	float bc_norm;
	float x,y,z;

	// The coordinate frame orign lies in the atom B for the ABCD dihedral torsion
	orign.x = p_atoms[*atmB-1].coord.x;
	orign.y = p_atoms[*atmB-1].coord.y;
	orign.z = p_atoms[*atmB-1].coord.z;

	bc.x = p_atoms[*atmC-1].coord.x - orign.x;
	bc.y = p_atoms[*atmC-1].coord.y - orign.y;
	bc.z = p_atoms[*atmC-1].coord.z - orign.z;

	bc_norm = sqrt( (bc.x*bc.x)+(bc.y*bc.y)+(bc.z*bc.z) );

	u.x = bc.x / bc_norm;
	u.y = bc.y / bc_norm;
	u.z = bc.z / bc_norm;

	x = p_atoms[*atmD-1].coord.x - orign.x;
	y = p_atoms[*atmD-1].coord.y - orign.y;
	z = p_atoms[*atmD-1].coord.z - orign.z;

	p_atoms[*atmD-1].coord.x =  ((u.x*u.x*(1-cos(*dih)))+( 1 *cos(*dih)))*x + ((u.x*u.y*(1-cos(*dih)))-(u.z*sin(*dih)))*y + ((u.x*u.z*(1-cos(*dih)))+(u.y*sin(*dih)))*z;
	p_atoms[*atmD-1].coord.y =  ((u.x*u.y*(1-cos(*dih)))+(u.z*sin(*dih)))*x + ((u.y*u.y*(1-cos(*dih)))+( 1 *cos(*dih)))*y + ((u.y*u.z*(1-cos(*dih)))-(u.x*sin(*dih)))*z;
	p_atoms[*atmD-1].coord.z =  ((u.x*u.z*(1-cos(*dih)))-(u.y*sin(*dih)))*x + ((u.z*u.y*(1-cos(*dih)))+(u.x*sin(*dih)))*y + ((u.z*u.z*(1-cos(*dih)))+( 1 *cos(*dih)))*z;

	p_atoms[*atmD-1].coord.x += orign.x;
	p_atoms[*atmD-1].coord.y += orign.y;
	p_atoms[*atmD-1].coord.z += orign.z;

}

/** Rotates all atoms of protein that are after num_res_first
* prot is the protein which will be moved (rotated)
* atomB fixed atom
* atomC fixed atom
* num_res_first number of residue that was rotated previously
*/
static void rotate_all_atoms(protein_t *prot, const int *atomB, const int *atomC, 
	const int *num_res_first, const float *angle){
	int forward_residue;
	forward_residue = *num_res_first + 1; // Obtain the number of next residue

	if (forward_residue <= prot->p_topol->numres){
		// From first atom of next_residue until last atom of protein, rotate all of them.
		for (int i = prot->p_topol->range_atoms[forward_residue-1].first_atom; i < prot->p_topol->numatom; i++){
			rotation_by_angle_dih(prot->p_atoms, atomB, atomC, &i, angle);
		}
	}

}


/** Rotates protein in a PSI dihedral movement
* prot is the protein which will be moved (rotated)
* num_res_first means the first residue will be moved. The dihedral moviment is continued in  
*               foward residue until last.
* angle is the value of rotated angle
*/
void rotation_psi(protein_t *prot, const int *num_res_first, const float *angle){

	//rotates first residue
	rotation_by_angle_dih(prot->p_atoms, 
		&prot->p_topol->psi[*num_res_first-1].fixed_atoms[0], 
		&prot->p_topol->psi[*num_res_first-1].fixed_atoms[1],
	    &prot->p_topol->psi[*num_res_first-1].moved_atoms[0], angle);

	// rotates all atoms from the forward residue
	rotate_all_atoms(prot, &prot->p_topol->psi[*num_res_first-1].fixed_atoms[0], 
		&prot->p_topol->psi[*num_res_first-1].fixed_atoms[1], 
		num_res_first, angle);
}


/** Rotates protein in a PHI dihedral movement
* prot is the protein which will be moved (rotated)
* num_res_first means the first residue will be moved. The dihedral moviment is continued in  
*               foward residue until last.
* angle is the value of rotated angle
*/
void rotation_phi(protein_t *prot, const int *num_res_first, const float *angle){

	//rotates all moved atoms of first residue
	for (int i = 0; i < prot->p_topol->phi[*num_res_first-1].num_moved; i++){
		rotation_by_angle_dih(prot->p_atoms, 
			&prot->p_topol->phi[*num_res_first-1].fixed_atoms[0], 
			&prot->p_topol->phi[*num_res_first-1].fixed_atoms[1],
	    	&prot->p_topol->phi[*num_res_first-1].moved_atoms[i], angle);
	}

	// rotates all atoms from the forward residue
	rotate_all_atoms(prot, &prot->p_topol->phi[*num_res_first-1].fixed_atoms[0], 
		&prot->p_topol->phi[*num_res_first-1].fixed_atoms[1], 
		num_res_first, angle);
}

/** Rotates protein in a OMEGA dihedral movement
* prot is the protein which will be moved (rotated)
* num_res_first means the first residue will be moved. The dihedral moviment is continued in  
*               foward residue until last.
* angle is the value of rotated angle
*/
void rotation_omega(protein_t *prot, const int *num_res_first, const float *angle){

	int num_res_next;

	//omega rotation is from forward residue of num_res_first
	num_res_next = *num_res_first + 1;

	/* Checking the number of residue. 
	* if num_res_next is greater than number of residue of protein 
	* no makes rotation in protein.
	*/
	if (num_res_next <= prot->p_topol->numres){
		//rotates all moved atoms of first residue
		for (int i = 0; i < prot->p_topol->omega[num_res_next-1].num_moved; i++){
			rotation_by_angle_dih(prot->p_atoms, 
				&prot->p_topol->omega[num_res_next-1].fixed_atoms[0], 
				&prot->p_topol->omega[num_res_next-1].fixed_atoms[1],
		    	&prot->p_topol->omega[num_res_next-1].moved_atoms[i], angle);
		}

		// rotates all atoms from the forward residue
		rotate_all_atoms(prot, &prot->p_topol->omega[num_res_next-1].fixed_atoms[0], 
			&prot->p_topol->omega[num_res_next-1].fixed_atoms[1], 
			&num_res_next, angle);
	}
}
