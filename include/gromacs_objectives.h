#ifndef OLD_GROMACS_OBJECTIVES_H
#define OLD_GROMACS_OBJECTIVES_H

#include "solution_types.h"
#include "parameters_type.h"

/** option_g_energy enum must be associated with type_fitness_energies.
 * This association is obtained througth
 * get_option_fitness_gromacs_from_type_fitness_energy function
 */
enum option_g_energy{gmx_potential_ener, gmx_edw_ener, gmx_elel_ener,
	gmx_hydrophobic, gmx_hydrophilic, gmx_total_area, gmx_gyrate,
	gmx_hbond, gmx_hbond_main, gmx_GBSA_Solvatation, gmx_GB_Polarization,
	gmx_Nonpolar_Sol, gmx_stride_total, gmx_stride_helix, gmx_stride_beta,
	gmx_ener_NR}; //gmx_hbond_side, gmx_hbond_side_main,
typedef enum option_g_energy option_g_energy_t;

typedef struct soption_fitness_gromacs{
	option_g_energy_t opt_g_energy;
	char value_opt[3];
	char option_name[30];
}option_fitness_gromacs_t;


void get_gromacs_objectives(solution_t *solutions, const input_parameters_t *in_para);
void init_gromacs_execution();
void finish_gromacs_execution();

#endif