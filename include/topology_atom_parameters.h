#include "topology_types.h"
#include "enums.h"


special_atom_parameters_t special_atom_parameters [] = {
		{atmCA_,atmCA,-1},
		{atmCA_plus,atmCA, +1},
        {atmC_,atmC,-1},
        {atmN_plus,atmN,+1},
        {atmO_,atmO,-1},
        {atmN_,atmN,-1}

};

topol_atoms_bond_parameters_t topol_atoms_bond_parameters [] = {
			{atmC,atmN,1.32},
			{atmC,atmN_plus,1.32},
			{atmCA,atmC,1.53},
			{atmCA,atmCB,1.53},
			{atmCA,atmN,1.45},
			{atmCA,atmHA,1.0},
			{atmCA,atmHA1,1.0},
			{atmCA,atmHA2,1.0},
			{atmC,atmO,1.43},
			{atmC,atmOT1,1.43},
			{atmC,atmOT2,1.43},
			{atmC,atmC,1.53},
			{atmC,atmHA,1.00},
			{atmC,atmHA1,1.00},
			{atmC,atmHA2,1.00},
			{atmCB,atmCA,1.53},
			{atmCB,atmCG,1.53},
			{atmNE,atmCD,1.46},
			{atmCZ,atmNE,1.329},
			{atmCZ,atmNH1,1.329},
			{atmNH2,atmCZ,1.326},
			{atmN,atmHN,1.0},
			{atmN,atmC_,1.329},
			{atmCB,atmHA,1.0},
			{atmCB,atmHB1,1.0},
			{atmCB,atmHB2,1.0},
			{atmCB,atmHB3,1.0},
			{atmCG,atmHG1,1.0},
			{atmCG,atmHG2,1.0},
			{atmCD,atmHD1,1.0},
			{atmCD,atmHD2,1.0},
			{atmCD,atmCG,1.53},
			{atmNE,atmHE,1.0},
			{atmNH1,atmHH11,1.0},
			{atmNH1,atmHH12,1.0},
			{atmNH2,atmHH21,1.0},
			{atmNH2,atmHH22,1.0},
			{atmN,atmH1, 1.0},
			{atmN,atmH2,1.0},
			{atmN,atmH3,1.0},
			{atmO,atmC,1.21}

};

/*The angle values are radius*/
topol_atoms_bond_angles_parameters_t topol_atoms_bond_angles_parameters [] = {
		{atmN,atmC_,atmCA_,2.0281}, //116.200
		{atmN,atmC_,atmCA,2.0281}, //116.200
		{atmN,atmCB,atmCA,2.0281}, //116.200
		{atmCA,atmN,atmC_,2.124}, //121.700
		{atmCB,atmCA,atmN,1.9269}, //110.400
		{atmCB,atmC,atmCA,1.9914}, //114.100
		{atmCG,atmCB,atmCA,1.9914}, //114.100
		{atmCD,atmCG,atmCB,1.9425},//111.300
		{atmC,atmCA,atmN,1.9409}, //111.200
		{atmC,atmCA,atmCB,1.9233}, //110.200
		{atmCZ,atmNE,atmCD,2.167}, // 124.200
		{atmNE,atmCD,atmCG,1.9545}, //112.000
		{atmNH2,atmCZ,atmNE,2.094}, //120.000
		{atmNH1,atmCZ,atmNE,2.094}, //120.000
		{atmNH1,atmCZ,atmNH2,2.094}, //120.000
		{atmO,atmC,atmN_plus,2.1468}, //123.0
		{atmO,atmC,atmCA,2.1083}, //120.800
		{atmOT1,atmC,atmCA,2.1083}, //120.800
		{atmOT1,atmC,atmN_plus,2.1468}, //123.0
		{atmOT2,atmC,atmOT1,3.1415}, //180.0
		{atmH1,atmCA,atmN,1.9111}, //109.50
		{atmH2,atmCA,atmN,1.9111}, //109.50
		{atmH3,atmCA,atmN,1.9111}, //109.50
		{atmHB1,atmCA,atmCB,1.9111}, //109.50
		{atmHB2,atmCA,atmCB,1.9111}, //109.50
		{atmHB3,atmCA,atmCB,1.9111}, //109.50
		{atmHA,atmN,atmCA,1.9111}, //109.50
		{atmHA1,atmN,atmCA,1.9111}, //109.50
		{atmHA2,atmN,atmCA,1.9111}, //109.50
		{atmHD1,atmCG,atmCD,1.9111}, //109.50
		{atmHD2,atmCG,atmCD,1.9111}, //109.50
		{atmHE,atmNE,atmCD,1.9111}, //109.50
		{atmHG1,atmCG,atmCB,1.9111}, //109.50
		{atmHG2,atmCG,atmCB,1.9111}, //109.50
		{atmHH11,atmNH1,atmCZ,1.9111}, //109.50
		{atmHH12,atmNH1,atmCZ,1.9111}, //109.50
		{atmHH21,atmNH2,atmCZ,1.9111}, //109.50
		{atmHH22,atmNH2,atmCZ,1.9111}, //109.50
		{atmHN,atmN,atmC,1.9111}, //109.50
		{atmHN,atmN,atmC_,1.9111}, //109.50
		{atmHN,atmN,atmCA,1.9111} //109.50
		//{atmOT2,atmC,atmCA,120.800},

};



/*
 * That parametres were based on:
 * Simon C. Lovell, J. Michael Word, Jane S. Richardson, and David C. Richardson,
 * The Penultimate Rotamer Library,
 * PROTEINS: Structure, Function, and Genetics 40:389 – 408 (2000)
 * */


