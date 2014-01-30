#ifndef OLD_TOPOLOGY_CHARMM27_PARAMETERS_H
#define OLD_TOPOLOGY_CHARMM27_PARAMETERS_H

#include"topology_types.h"
#include"aminoacids/charmm27/GLY.h"
#include"aminoacids/charmm27/ALA.h"
#include"aminoacids/charmm27/ARG.h"
#include"aminoacids/charmm27/ASN.h"
#include"aminoacids/charmm27/ASP.h"
#include"aminoacids/charmm27/CYS.h"
#include"aminoacids/charmm27/GLN.h"
#include"aminoacids/charmm27/GLU.h"
#include"aminoacids/charmm27/ILE.h"
#include"aminoacids/charmm27/LEU.h"
#include"aminoacids/charmm27/LYS.h"
#include"aminoacids/charmm27/SER.h"
#include"aminoacids/charmm27/MET.h"
#include"aminoacids/charmm27/PHE.h"
#include"aminoacids/charmm27/PRO.h"
#include"aminoacids/charmm27/THR.h"
#include"aminoacids/charmm27/TRP.h"
#include"aminoacids/charmm27/TYR.h"
#include"aminoacids/charmm27/VAL.h"
#include"aminoacids/charmm27/HIS.h"
#include"aminoacids/charmm27/HSE.h"
#include"aminoacids/charmm27/HSD.h"

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_omega [] = {
		{atmCA, atmC, atmN_plus, atmCA_plus} //atmHA
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_psi [] = {
		{atmN, atmCA, atmC, atmN_plus} //atmHA
};

topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_phi [] = {
		{atmC_, atmN, atmCA, atmC}
};

/*I modified my phi C-Terminal. I am not sure if I done correctly modification
 * Basically, I changed atmN_plus to atmOT1, because atmN_plus does not exists.*/
topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_phi_C_Terminal [] = {
		{atmC_, atmN, atmCA, atmC}
};

//topol_residue_atoms_dihedral_angles_t residue_atom_diehdral_phi_C_Terminal [] = {
		//{atmC, atmCA, atmN, atmOT1}//Oxygen is not present in C-Term. It was changed by OT1
//};
/*
 * ATENCAO: ALA foi marcado sendo 1 cadeia lateral
 *         ARG segundo a biblioteca de rotameros do gui tem 5. Mas o protpred
 *         marca 4
*/

/*The index of topol_residues_ff is type_aminos_t when N-Terminal*/
topol_residues_t topol_residues_ff_N_Terminal []= {
		                               {aALA,"ALA","A",12,1,1,residue_atoms_ALA_N_Terminal,residue_atoms_bonds_ALA_N_Terminal,5,residue_atoms_bonds_angles_ALA_N_Terminal,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi, residue_atom_diehdral_rot_ALA, 2,residue_atoms_dihedral_angles_type_ALA_N_Terminal, residue_atom_diehdral_omega},
		                               {aARG,"ARG","R",26,4,0,residue_atoms_ARG_N_Terminal,residue_atoms_bonds_ARG_N_Terminal,9,residue_atoms_bonds_angles_ARG,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi, residue_atom_diehdral_rot_ARG, 3, residue_atoms_dihedral_angles_type_ARG, residue_atom_diehdral_omega},
		                               {aGLY,"GLY","G",9,0,1,residue_atoms_GLY_N_Terminal,residue_atoms_bonds_GLY_N_Terminal,4,residue_atoms_bonds_angles_GLY,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi,NULL,2,residue_atoms_dihedral_angles_type_GLY, residue_atom_diehdral_omega},
		                               {aASN,"ASN","N",16,2,0,residue_atoms_ASN_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ASN,2,NULL, residue_atom_diehdral_omega},
		                               {aASP,"ASP","D",14,2,0,residue_atoms_ASP_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ASP,2,NULL, residue_atom_diehdral_omega},
		                               {aCYS,"CYS","C",13,1,0,residue_atoms_CYS_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_CYS,2,NULL, residue_atom_diehdral_omega},
		                               {aGLN,"GLN","Q",19,3,0,residue_atoms_GLN_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_GLN,2,NULL, residue_atom_diehdral_omega},
		                               {aGLU,"GLU","E",17,3,0,residue_atoms_GLU_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_GLU,2,NULL, residue_atom_diehdral_omega},
		                               {aILE,"ILE","I",21,2,0,residue_atoms_ILE_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ILE,2,NULL, residue_atom_diehdral_omega},
		                               {aLEU,"LEU","L",21,2,0,residue_atoms_LEU_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_LEU,2,NULL, residue_atom_diehdral_omega},
		                               {aLYS,"LYS","K",24,4,0,residue_atoms_LYS_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_LYS,2,NULL, residue_atom_diehdral_omega},
		                               {aSER,"SER","S",13,1,0,residue_atoms_SER_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_SER,2,NULL, residue_atom_diehdral_omega},
		                               {aMET,"MET","M",19,3,0,residue_atoms_MET_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_MET,2,NULL, residue_atom_diehdral_omega},
		                               {aPHE,"PHE","F",22,2,0,residue_atoms_PHE_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_PHE,2,NULL, residue_atom_diehdral_omega},
		                               {aPRO,"PRO","P",16,2,0,residue_atoms_PRO_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_PRO,2,NULL, residue_atom_diehdral_omega},
		                               {aTHR,"THR","T",16,1,0,residue_atoms_THR_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_THR,2,NULL, residue_atom_diehdral_omega},
		                               {aTRP,"TRP","W",26,2,0,residue_atoms_TRP_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_TRP,2,NULL, residue_atom_diehdral_omega},
		                               {aTYR,"TYR","Y",23,2,0,residue_atoms_TYR_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_TYR,2,NULL, residue_atom_diehdral_omega},
		                               {aVAL,"VAL","V",18,1,0,residue_atoms_VAL_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_VAL,2,NULL, residue_atom_diehdral_omega},
		                               {aHIS,"HIS","H",19,2,0,residue_atoms_HIS_N_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_HIS,2,NULL, residue_atom_diehdral_omega}


};

/*The index of topol_residues_ff is type_aminos_t */
topol_residues_t topol_residues_ff []= {
		                               {aALA,"ALA","A",10,1,1,residue_atoms_ALA,residue_atoms_bonds_ALA,7,residue_atoms_bonds_angles_ALA,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ALA, 2,residue_atoms_dihedral_angles_type_ALA, residue_atom_diehdral_omega},
		                               {aARG,"ARG","R",24,4,0,residue_atoms_ARG,residue_atoms_bonds_ARG,11,residue_atoms_bonds_angles_ARG,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi, residue_atom_diehdral_rot_ARG, 3, residue_atoms_dihedral_angles_type_ARG, residue_atom_diehdral_omega},
		                               {aGLY,"GLY","G",7,0,1,residue_atoms_GLY,residue_atoms_bonds_GLY,4,residue_atoms_bonds_angles_GLY,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi,NULL,2,residue_atoms_dihedral_angles_type_GLY, residue_atom_diehdral_omega},
		                               {aASN,"ASN","N",14,2,0,residue_atoms_ASN,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ASN,2,NULL, residue_atom_diehdral_omega},
		                               {aASP,"ASP","D",12,2,0,residue_atoms_ASP,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ASP,2,NULL, residue_atom_diehdral_omega},
		                               {aCYS,"CYS","C",10,1,0,residue_atoms_CYS,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_CYS,2,NULL, residue_atom_diehdral_omega},
		                               {aGLN,"GLN","Q",17,3,0,residue_atoms_GLN,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_GLN,2,NULL, residue_atom_diehdral_omega},
		                               {aGLU,"GLU","E",15,3,0,residue_atoms_GLU,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_GLU,2,NULL, residue_atom_diehdral_omega},
		                               {aILE,"ILE","I",19,2,0,residue_atoms_ILE,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ILE,2,NULL, residue_atom_diehdral_omega},
		                               {aLEU,"LEU","L",19,2,0,residue_atoms_LEU,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_LEU,2,NULL, residue_atom_diehdral_omega},
		                               {aLYS,"LYS","K",22,4,0,residue_atoms_LYS,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_LYS,2,NULL, residue_atom_diehdral_omega},
		                               {aSER,"SER","S",11,1,0,residue_atoms_SER,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_SER,2,NULL, residue_atom_diehdral_omega},
		                               {aMET,"MET","M",17,3,0,residue_atoms_MET,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_MET,2,NULL, residue_atom_diehdral_omega},
		                               {aPHE,"PHE","F",20,2,0,residue_atoms_PHE,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_PHE,2,NULL, residue_atom_diehdral_omega},
		                               {aPRO,"PRO","P",14,2,0,residue_atoms_PRO,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_PRO,2,NULL, residue_atom_diehdral_omega},
		                               {aTHR,"THR","T",14,1,0,residue_atoms_THR,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_THR,2,NULL, residue_atom_diehdral_omega},
		                               {aTRP,"TRP","W",24,2,0,residue_atoms_TRP,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_TRP,2,NULL, residue_atom_diehdral_omega},
		                               {aTYR,"TYR","Y",21,2,0,residue_atoms_TYR,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_TYR,2,NULL, residue_atom_diehdral_omega},
		                               {aVAL,"VAL","V",16,1,0,residue_atoms_VAL,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_VAL,2,NULL, residue_atom_diehdral_omega},
		                               {aHIS,"HIS","H",17,2,0,residue_atoms_HIS,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_HIS,2,NULL, residue_atom_diehdral_omega},
		                               {aHSE,"HIS","H",17,2,0,residue_atoms_HSE,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_HSE,2,NULL, residue_atom_diehdral_omega},
		                               {aHSD,"HIS","H",17,2,0,residue_atoms_HSD,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_HSD,2,NULL, residue_atom_diehdral_omega}


};

/*The index of topol_residues_ff is type_aminos_t when C-Terminal*/
topol_residues_t topol_residues_ff_C_Terminal []= {
        {aALA,"ALA","A",11,1,1,residue_atoms_ALA_C_Terminal,residue_atoms_bonds_ALA_C_Terminal,5,residue_atoms_bonds_angles_ALA_N_Terminal,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi, residue_atom_diehdral_rot_ALA, 2,residue_atoms_dihedral_angles_type_ALA_C_Terminal, residue_atom_diehdral_omega},
		{aARG,"ARG","R",25,4,0,residue_atoms_ARG_C_Terminal,residue_atoms_bonds_ARG_C_Terminal,11,residue_atoms_bonds_angles_ARG_C_Terminal,4, residue_atom_diehdral_phi_C_Terminal, residue_atom_diehdral_psi, residue_atom_diehdral_rot_ARG, 3,residue_atoms_dihedral_angles_type_ARG_C_Terminal, residue_atom_diehdral_omega},
		{aGLY,"GLY","G",8,0,1,residue_atoms_GLY_C_Terminal,residue_atoms_bonds_GLY_C_Terminal,4,residue_atoms_bonds_angles_GLY,4, residue_atom_diehdral_phi, residue_atom_diehdral_psi,NULL,2,residue_atoms_dihedral_angles_type_GLY_C_Terminal, residue_atom_diehdral_omega},
		{aASN,"ASN","N",15,2,0,residue_atoms_ASN_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ASN,2,NULL, residue_atom_diehdral_omega},
		{aASP,"ASP","D",13,2,0,residue_atoms_ASP_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ASP,2,NULL, residue_atom_diehdral_omega},
		{aCYS,"CYS","C",12,1,0,residue_atoms_CYS_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_CYS,2,NULL, residue_atom_diehdral_omega},
		{aGLN,"GLN","Q",18,3,0,residue_atoms_GLN_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_GLN,2,NULL, residue_atom_diehdral_omega},
		{aGLU,"GLU","E",16,3,0,residue_atoms_GLU_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_GLU,2,NULL, residue_atom_diehdral_omega},
        {aILE,"ILE","I",20,2,0,residue_atoms_ILE_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_ILE,2,NULL, residue_atom_diehdral_omega},
        {aLEU,"LEU","L",20,2,0,residue_atoms_LEU_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_LEU,2,NULL, residue_atom_diehdral_omega},
        {aLYS,"LYS","K",23,4,0,residue_atoms_LYS_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_LYS,2,NULL, residue_atom_diehdral_omega},
        {aSER,"SER","S",12,1,0,residue_atoms_SER_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_SER,2,NULL, residue_atom_diehdral_omega},
        {aMET,"MET","M",18,3,0,residue_atoms_MET_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_MET,2,NULL, residue_atom_diehdral_omega},
        {aPHE,"PHE","F",21,2,0,residue_atoms_PHE_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_PHE,2,NULL, residue_atom_diehdral_omega},
        {aPRO,"PRO","P",15,2,0,residue_atoms_PRO_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_PRO,2,NULL, residue_atom_diehdral_omega},
        {aTHR,"THR","T",15,1,0,residue_atoms_THR_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_THR,2,NULL, residue_atom_diehdral_omega},
        {aTRP,"TRP","W",25,2,0,residue_atoms_TRP_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_TRP,2,NULL, residue_atom_diehdral_omega},
        {aTYR,"TYR","Y",22,2,0,residue_atoms_TYR_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_TYR,2,NULL, residue_atom_diehdral_omega},
        {aVAL,"VAL","V",17,1,0,residue_atoms_VAL_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_VAL,2,NULL, residue_atom_diehdral_omega},
        {aHIS,"HIS","H",18,2,0,residue_atoms_HIS_C_Terminal,NULL,8,NULL,4,residue_atom_diehdral_phi, residue_atom_diehdral_psi,residue_atom_diehdral_rot_HIS,2,NULL, residue_atom_diehdral_omega}
};
#endif
