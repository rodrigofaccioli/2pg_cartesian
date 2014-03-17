#ifndef OLD_TOPOLOGY_TYPES_H
#define OLD_TOPOLOGY_TYPES_H

/** represents the resideue and atoms
* num_fixed number of fixed atoms. Atoms that will be fixed.
* fixed_atoms array of index fixed atoms.
* num_moved number of moved atoms. Atoms that will be moved.
* moved_atoms array of index moved atoms.
*/
typedef struct stop_residue_atom_info{
	int num_fixed;
	int *fixed_atoms;
	int num_moved;
	int *moved_atoms;
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


/** represents the information about side chain of residue
* num_chi means the number of chi that the residue has. This number is assigned
* by get_number_chi function at topology.c file.
* atoms_chi is an array of atoms (fixed and moved) for each chi. In other words,
* the elements number of this array is given by num_chi.
*/
typedef struct stop_residue_side_chains{
	int num_chi;
	top_residue_atom_info_t *atoms_chi;
}top_residue_side_chains_t;


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
	top_residue_atom_info_t *omega;
	top_residue_side_chains_t *side_chains;
 }top_global_t;



#endif
