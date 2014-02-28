#ifndef OLD_PARAMETERS_TYPE_H
#define OLD_PARAMETERS_TYPE_H

#include "enums.h"

typedef struct sinput_parameters_extended_chain{
	char *local_execute;
	char *seq_protein_file_name;
	int res_reference;
	char *file_final_pdb;
	char *path_database;
	int steps;
}input_parameters_extended_chain_t;

typedef struct sinput_parameters
{
    int size_population;
    char *seq_protein_file_name;
    char *top_file;
    char *z_matrix_file;
    char *initial_pop_file_name;
    int number_fitness;
    int number_generation;
    int number_individual_select_reproduce;
    float crossover_rate;
    float mutation_rate;
    float blx_cros_rate;
    float point_1_cros_rate;
    float point_2_cros_rate;
    int number_crossover; //it shows how many crossovers the user chose. This parameter is obtained by set_parameter_number_crossover function
    /*It stores the crossovers the user chose. It is considered a crossover chose when its rate is greater than 0.
     * Otherwise, the crossover is not chosen. For example, point_1_cros_rate > 0, in crossovers contains at least crossoer_point_1.
     */
    type_crossoers_t *crossovers;
    type_energy_minimization_t gromacs_energy_min; //Set which energy minimization: none, implicit or explicit. none is default.
    type_energy_minimization_t gromacs_energy_min_gen_oper; ////Set which energy minimization for genetic operators: none, implicit or explicit. none is default.
    int processor_number; //Set the number of processor to run. Default is 1. Otherwise, it will be considered run in parallel
    char *file_final_pdb;
    char *path_file_native_protein;
    char *path_database;
    char *path_program_compute_energy; //script to run simulation. Ex: run_gromacs_compute_energy_potential.sh
    char *path_program_minimization_program; //script to minimization energy. Calls at build_random_protein. Ex: run_gromacs_minimization.py
    char *path_program_g_energy; //script calls to run g_energy program. Ex: run_g_energy.sh
    char *path_program_clean_simulation; //script calls clean simulation. Ex: clean_simulation.sh
    char *path_program_HIS_protonation; //script calls HIS protonation. Ex: run_his_protonation.sh
    char *path_program_rmsd; //script calls rmsd calculation. Ex: run_g_rms.sh
    char *path_local_execute;
    char *path_gromacs_programs;
    char *computed_energies_gromacs_file;
    char *energy_file_xvg;
    char *path_program_read_energy;
    char *path_program_read_xvg_rmsd;
    char *path_program_g_sas;
    char *path_program_read_g_sas;
    char *path_program_g_gyrate;
    char *path_program_read_g_gyrate;
    char *path_program_g_hbond;
    char *path_program_read_g_hbond;
    char *computed_energy_value_file;
    char *computed_areas_g_sas_file;
    char *computed_radius_g_gyrate_file;
    char *computed_g_hbond_file;
    char *path_program_stride;
    /* This vector stores the kind of energies which will be used to obtain
     * the fitness value.
     * The number of elements must be equal number_fitness parameter. Otherwise,
     * a fatal error is executed.
     */
    type_fitness_energies_t *fitness_energies;
    /* This vector stores the weights for the weighted average. It is used at
     * MEAT algorithm for example.
     * The number of elements must be equal number_fitness parameter. Otherwise,
     * a fatal error is executed.
     * It is important to know its parameter works together with fitness_energies
     * parameter. If the number of values at fitness_energies are 3, this value
     * will be equal at weights_fitness.
     */
    float *weights_fitness;
    float blx_alfa;
    int number_archive; //Indicates the number of archives that SPEA2 algorithm works
    type_rotamer_library_t rotamer_library; // shows what kind of rotamer library
    type_objective_analysis_t objective_analysis; // shows what kind of objective analysis
    char *path_dimo_sources; // shows where DIMO source is
    char *path_call_GreedyTreeGenerator2PG;
    float max_mutation_range;
    float individual_mutation_rate;
//    float rate_mutation_phi;
    float min_angle_mutation_phi;
    float max_angle_mutation_phi;
//    float rate_mutation_psi;
    float min_angle_mutation_psi;
    float max_angle_mutation_psi;
//    float rate_mutation_omega;
    float min_angle_mutation_omega;
    float max_angle_mutation_omega;

 }input_parameters_t;

#endif
