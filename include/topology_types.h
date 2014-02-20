#ifndef OLD_TOPOLOGY_TYPES_H
#define OLD_TOPOLOGY_TYPES_H

/** represents the resideue and atoms
* num_fixed number of fixed atoms. Atoms that will be fixed.
* fixed_atoms array of index fixed atoms.
* num_moved number of moved atoms. Atoms that will be moved.
* moved_atoms array of index moved atoms.
* num_side_chains number of atoms in side chains.
* atoms_side_chains array of index atoms in side chains.
*/
typedef struct stop_residue_atom_info{
	int num_fixed;
	int *fixed_atoms;
	int num_moved;
	int *moved_atoms;
	int num_side_chains;
	int *atoms_side_chains;
}top_residue_atom_info_t;

/** represents the range of atoms
* first_atom number of first atom of residue
* last_atom number of last atom of residue
* Important: number means the value of PDB. This number is not started with 0.
* Therefore, if necessary access the index of atoms, it have to decrement 1 like 
* first_atom -1.
*/
typedef struct stop_residue_range_atoms{
	int first_atom;
	int last_atom;
}top_residue_range_atoms_t;


/** represents all information about topology
* phi information about atoms for phi.
* psi information about atoms for psi.
* omega information about atoms for omega. It is not necessary since
*       the atoms will be moved (rotated) are all atoms from next residue.
* side_chains  information about atoms for side_chains.
*/
typedef struct stop_global{
	int numatom;
	int numres;
	top_residue_range_atoms_t *range_atoms;
	top_residue_atom_info_t *phi;
	top_residue_atom_info_t *psi;
	//top_residue_atom_info_t *omega;
	top_residue_atom_info_t *side_chains;
 }top_global_t;



#endif
