#ifndef OLD_ENUM_H
#define OLD_ENUM_H

//enum type_aminos {aGLY, aALA, aVAL, aLEU, aILE, aSER, aTHR, aCYS, aCYX, aPRO, aPHE, aTYR, aTRP, aHIS, aHID, aHIE, aASP, aASN, aGLU, aGLN,
	              //aMET, aLYS, aARG, aORN, aAIB, aPCA, aFOR, aACE, aNH2, aNME, aUNK, aNR};

enum boolean {bfalse, btrue};
typedef enum boolean boolean_t;

enum type_aminos {aGLY, aALA, aARG, aVAL, aLEU, aILE, aSER, aTHR, aCYS, aPRO, aPHE, aTYR, aHIS, aASP, aASN, aGLU, aMET, aLYS,
	               aGLN, aTRP, aHSE, aHSD,
	               aNR /*This line must be the last*/
                 };
typedef enum type_aminos type_aminos_t;

/*All atoms which are supported by topology. These index are
 * applied to build the system topology
 *
 * 1) The atoms tmC_, atmCA_, atmCA_plus and atmN_plus
 * are special atoms. These atoms show that they are
 * in different residue. However, they represents the same
 * kind of atom. They are used to build the topology of
 * angle bond and bond length.
 * _ represents before residue
 * _plus represents next residue
 *
 */
enum type_atoms { atmC, atmCA, atmCB, atmCG,atmCG1, atmCG2,
	             atmCD, atmCD1, atmCD2,
	             atmCE, atmCE1, atmCE2, atmCE3,
	             atmCH2,
	             atmCZ, atmCZ2, atmCZ3,
	             atmN, atmND1, atmND2, atmNH1, atmNH2,
	             atmNE, atmNE1, atmNE2,
	             atmNZ,
	             atmO, atmOD1, atmOD2, atmOE1, atmOE2, atmOH,
	             atmOG, atmOG1,
	             atmH, atmHA, atmHA1, atmHA2, atmHA3,
	             atmHE, atmHE1,atmHE2, atmHE3, atmHE21, atmHE22,
	             atmHG, atmHG1, atmHG2, atmHG3,atmHG11, atmHG12, atmHG13,
	             atmHG21, atmHG22, atmHG23,
	             atmHN,
	             atmH1, atmH2, atmH3,
	             atmHB, atmHB1,atmHB2,atmHB3,
	             atmHD1, atmHD2, atmHD3, atmHD11, atmHD12, atmHD13, atmHD21,
	             atmHD22, atmHD23,
	             atmHH, atmHH11, atmHH12, atmHH2, atmHH21, atmHH22,
	             atmHZ, atmHZ1, atmHZ2, atmHZ3,
	             atmC_, atmCA_, atmCA_plus, atmN_plus, atmO_, atmN_,
	             atmOT1, atmOT2, atmOXT,
	             atmSG,atmSD,
	             atmNR /*This line must be the last*/
                 };
typedef enum type_atoms type_atoms_t;

enum default_messages {msgOPTION_WRONG, msgNR};
typedef enum default_messages default_messages_t;

/* Shows the kind of dihedral angles which are supported by algorithm.
 * It is used in nerf topology for example
 * angl_typ_dieh_180 means the value of dihedral angle is 180 degrees.
 * angl_typ_dieh_0 represents the dihedral angle makes by CD-NE-CZ-NH1
 * angl_typ_dieh_1 represents the dihedral angle makes by N-C-CA-O
 * angl_typ_dieh_2 represents the dihedral angle makes by N-CA-C-CB
 * angl_typ_dieh_3 represents the dihedral angle makes by randomly
 * angl_typ_dieh_4 represents chi2 - 180 or chi2 + 180. Depends on value of chi2
 * angl_typ_dieh_5 represents chi3 - 180 or chi3 + 180. Depends on value of chi3
 * angl_typ_dieh_6 represents chi1 - 120
 * angl_typ_dieh_7 represents chi2 - 120
 * angl_typ_dieh_trans_123_ represents -123.0 degrees
 * dihedral (RD): RD - 180 if RD is positive. Otherwise, 180 + RD
 */
enum type_dihedral_angles {angl_phi, angl_psi,angl_psi_,
	angl_chi1, angl_chi2, angl_chi3, angl_chi4, angl_typ_dieh_1,
	angl_typ_dieh_2, angl_typ_dieh_3,
	angl_typ_dieh_4, angl_typ_dieh_5, angl_typ_dieh_6, angl_typ_dieh_7,
	angl_typ_dieh_180, angl_typ_dieh_omega, angl_typ_dieh_0, angl_typ_dieh_90,
	angl_typ_dieh_117_, angl_typ_dieh_trans_120, angl_typ_dieh_trans_123_,
	angl_typ_dieh_trans_240, angl_typ_dieh_29_6, angl_typ_dieh_37_4_};
typedef enum type_dihedral_angles type_dihedral_angles_t;

/*Represents the energies which are used to obtain the fitness.*/
enum type_fitness_energies{fit_ener_potential, fit_ener_edw, fit_ener_ele,
	fit_hydrophobic, fit_hydrophilic, fit_total_area, fit_gyrate,
	fit_hbond, fit_hbond_main, fit_GBSA_Solvatation, fit_stride_total,
	fit_stride_helix, fit_stride_beta, fit_ener_NR}; //fit_hbond_side, fit_hbond_side_main,
typedef enum type_fitness_energies type_fitness_energies_t;

/*Represents kind of energy minimization is used.*/
enum type_energy_minimization{ener_min_none, ener_min_implicit, ener_min_explicit,
	ener_min__NR};
typedef enum type_energy_minimization type_energy_minimization_t;

/*Represents kind of rotamer library.*/
enum type_rotamer_library{rotamer_library_none, rotamer_library_cad_tuffery,
	rotamer_library_NR};
typedef enum type_rotamer_library type_rotamer_library_t;

/*Represents kind of crossover */
enum type_crossoers{crossoer_point_1, crossoer_point_2, crossoer_blx_alpha,
	crossoer_NR};
typedef enum type_crossoers type_crossoers_t;

/*Represents kind of objective_analysis */
enum type_objective_analysis {objective_analysis_none,
	objective_analysis_dimo_GreedyTreeGenerator, objective_analysis_NR};

typedef enum type_objective_analysis type_objective_analysis_t;

#endif





