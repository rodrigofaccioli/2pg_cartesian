#ifndef OLD_GROMACS_OBJECTIVES_H
#define OLD_GROMACS_OBJECTIVES_H

#include "solution_types.h"
#include "parameters_type.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

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

_2PG_CARTESIAN_EXPORT
void get_gromacs_objectives_of_solution(solution_t *solution, 
	const input_parameters_t *in_para, const int *ind);
_2PG_CARTESIAN_EXPORT
void get_gromacs_objectives(solution_t *solutions, const input_parameters_t *in_para);
_2PG_CARTESIAN_EXPORT
void init_gromacs_execution();
_2PG_CARTESIAN_EXPORT
void finish_gromacs_execution();
void build_tpr_file(const char *pdbfile, const char *local_execute,
		const char *path_gromacs_programs, const char *force_field, const char *mdp_file);
void call_pdb2gmx_for_pattern_atom_names(const char *pdbfile, const char *local_execute,
		const char *path_gromacs_programs, const char *force_field);
void clean_gromacs_simulation(const char *path_local_execute);
#endif
