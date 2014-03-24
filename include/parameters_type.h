#ifndef OLD_PARAMETERS_TYPE_H
#define OLD_PARAMETERS_TYPE_H

#include "enums.h"


typedef struct sinput_parameters{
    char *path_local_execute;
    char *path_gromacs_programs;
    int size_population;
    char *seq_protein_file_name;    
    char *initial_pop_file_name;
    int number_fitness;
    type_fitness_energies_t *fitness_energies;
    int number_generation;
    float point_1_cros_rate;
    float point_2_cros_rate;
    int number_crossover; //it shows how many crossovers the user chose. This parameter is obtained by set_parameter_number_crossover function
    /*It stores the crossovers the user chose. It is considered a crossover chose when its rate is greater than 0.
     * Otherwise, the crossover is not chosen. For example, point_1_cros_rate > 0, in crossovers contains at least crossoer_point_1.
     */
    type_crossoers_t *crossovers;
    type_energy_minimization_t gromacs_energy_min; //Set which energy minimization: none, implicit or explicit. none is default.
    type_energy_minimization_t gromacs_energy_min_gen_oper; ////Set which energy minimization for genetic operators: none, implicit or explicit. none is default.
    char *computed_energies_gromacs_file;
    char *energy_file_xvg;
    char *computed_energy_value_file;
    char *computed_areas_g_sas_file;
    char *computed_radius_g_gyrate_file;
    char *computed_g_hbond_file;
    type_objective_analysis_t objective_analysis; // shows what kind of objective analysis
    char *path_dimo_sources; // shows where DIMO source is
    char *path_call_GreedyTreeGenerator2PG;
    float individual_mutation_rate;
    float min_angle_mutation_phi;
    float max_angle_mutation_phi;
    float min_angle_mutation_psi;
    float max_angle_mutation_psi;
    float min_angle_mutation_omega;
    float max_angle_mutation_omega;
    float min_angle_mutation_side_chain;
    float max_angle_mutation_side_chain;
    /*Informs the force field*/
    char *force_field;
    /*Informs the mdp file name*/
    char *mdp_file;
    /* Number of steps to Monte Carlo Algorithm */
    int MonteCarloSteps;
 }input_parameters_t;

#endif
