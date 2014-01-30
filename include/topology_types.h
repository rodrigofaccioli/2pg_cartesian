#ifndef OLD_TOPOLOGY_TYPES_H
#define OLD_TOPOLOGY_TYPES_H

#include "enums.h"


/*Represent the atoms which are different residue.
 * They are used to obtain the bond length and
 * bond angle values.
 *  */
typedef struct sspecial_atom_parameters{
	type_atoms_t sp_atom_id;
    type_atoms_t atom_id;
    int residue_position;
}special_atom_parameters_t;

typedef struct stopol_atoms_bond_parameters{
	type_atoms_t atom_id_1;
	type_atoms_t atom_id_2;
	float bond_length;
} topol_atoms_bond_parameters_t;

typedef struct stopol_atoms_bond_angles_parameters{
	type_atoms_t atom_id_1;
	type_atoms_t atom_id_2;
	type_atoms_t atom_id_3;
	float angle_value;
} topol_atoms_bond_angles_parameters_t;


/*Represents how atoms of residue are bonds*/
typedef struct stopol_residue_atoms_bonds{
	type_atoms_t atom_id_1;
	type_atoms_t atom_id_2;
} topol_residue_atoms_bonds_t;

/*Represents the atoms of residue make angles*/
typedef struct stopol_residue_atoms_bonds_angles{
	type_atoms_t atom_id_1;
	type_atoms_t atom_id_2;
	type_atoms_t atom_id_3;
} topol_residue_atoms_bonds_angles_t;


/*Represents the atoms of residue which make
 * a dihedral angle. Dihedral angle means a angle formed by four atoms.
 */
typedef struct stopol_residue_atoms_dihedral_angles{
	type_atoms_t atom_id_1;
	type_atoms_t atom_id_2;
	type_atoms_t atom_id_3;
	type_atoms_t atom_id_4;
} topol_residue_atoms_dihedral_angles_t;

/*Represents the atoms of residue which make  a dihedral angle and its type.
 * Dihedral angle means an angle formed by four atoms. Its type is defined at
 * type_dihedral_angles structure.
 */
typedef struct stopol_residue_atoms_dihedral_angles_type{
	type_atoms_t atom_id_1;
	type_atoms_t atom_id_2;
	type_atoms_t atom_id_3;
	type_atoms_t atom_id_4;
	type_dihedral_angles_t type_dihedral;
} topol_residue_atoms_dihedral_angles_type_t;

typedef struct stopol_residue_atoms{
	type_atoms_t atom_id;
	const char *atom_name;
	const char *type;
	float charge;
} topol_residue_atoms_t;


typedef struct stopol_residues{
	const type_aminos_t res_id;
	const char *res_name;
	const char *res_name_1;
	const int nr_atoms;
	const int nr_side_chains;
	const int HP;
	const topol_residue_atoms_t *residue_atoms;
	const topol_residue_atoms_bonds_t *residue_atoms_bond;
	const int number_bond_angle;
	const topol_residue_atoms_bonds_angles_t *residue_atoms_bonds_angles;
	const int number_dihedral_angle;
	const topol_residue_atoms_dihedral_angles_t *residue_atoms_phi;
	const topol_residue_atoms_dihedral_angles_t *residue_atoms_psi;	
	const topol_residue_atoms_dihedral_angles_t *residue_atoms_side_chains;
	const int number_dihedral_angle_type;
	const topol_residue_atoms_dihedral_angles_type_t *residue_atoms_dihedral_angles_type;	
	const topol_residue_atoms_dihedral_angles_t *residue_atoms_omega;
} topol_residues_t;


typedef struct stop_global_atom{
	type_atoms_t atom_id;
	int atom_number;
	char atom_name[5];
	int res_number;
	type_aminos_t amino_id;
	char res_name[4];
	float charge;
} top_global_atom_t;

typedef struct stop_global_dihedral{
	int res_number;
	int atom_number1;
	int atom_number2;
	int atom_number3;
	int atom_number4;
} top_global_dihedral_t;

/*Represents the dihedral angles of protein and its type */
typedef struct stop_global_dihedral_angles_type{
	int res_number;
	int atom_number1;
	int atom_number2;
	int atom_number3;
	int atom_number4;
	type_dihedral_angles_t type_dihedral;
} top_global_dihedral_angles_type_t;

typedef struct stop_global_dihedral_side_chain{
	int res_number;
	int chi; //chi from 1 until 5
	int atom_number1;
	int atom_number2;
	int atom_number3;
	int atom_number4;
} top_global_dihedral_side_chain_t;

/*For each residue informs first and last atoms*/
typedef struct stop_global_res_atm{
	int res_number;
	int atom_first;
	int atom_last;
}top_global_res_atm_t;

typedef struct stop_global_res_atms_bond{
	int res_number;
	int atom_number1;
	int atom_number2;
	float bond_value;
}top_global_res_atms_bond_t;


typedef struct stop_global_res_atms_bond_angle{
	int res_number;
	int atom_number1;
	int atom_number2;
	int atom_number3;
	float angle_value;
}top_global_res_atms_bond_angle_t;


/*Represents the backbone of protein*/
typedef struct sprotein_backbone{
	int res_number;
	int atom_C;
	int atom_Ca;
	int atom_N;
	int atom_C_;
	int atom_N_plus;
}protein_backbone_t;


typedef struct stop_global{
	int numatom;
	int numres;
	int number_bond_angles;
	int number_protein_side_chains;
	int number_dihedral_angles_type;
	float protein_charge;
	top_global_atom_t *top_global_atom;
	top_global_res_atm_t *top_global_res_atm;
	top_global_res_atms_bond_t *top_global_res_atms_bond;
	top_global_res_atms_bond_angle_t *top_global_res_atms_bond_angle;
	top_global_dihedral_t *top_global_dieh_phi;
	top_global_dihedral_t *top_global_dieh_psi;
	top_global_dihedral_side_chain_t *top_global_dieh_side_chains;
	top_global_dihedral_angles_type_t *top_global_dihedral_angles_type;
	top_global_dihedral_t *top_global_dieh_omega;
}top_global_t;


#endif
