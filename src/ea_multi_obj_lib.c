#include <stdlib.h>

#include "ea_multi_obj_lib.h"
#include "energy_gromacs.h"
#include "defines.h"
#include "enums.h"
#include "messages.h"
#include "string_owner.h"
#include "ea.h"
#include "pdbio.h"
#include "functions.h"

void _fitness_gromacs(protein **population_p,
		const input_parameters_t *in_para, const top_global_t *top){
	int i;
	pdb_atom_t *pdb_atoms=NULL;
	char pdbfile_aux[20];
	char aux [10];
	char msg [50];
	int charge_value;
	char c_charge_value[15];
	option_g_energy *opt_fitness=NULL;

  pdb_atoms = allocate_pdbatom(&top->numatom);

	opt_fitness = Malloc(option_g_energy,in_para->number_fitness);
	for (i = 0; i < in_para->number_fitness; i++){
		opt_fitness[i] = get_option_g_energy_t_from_type_fitness_energy(&in_para->fitness_energies[i]);
	}
	charge_value = get_charge_topol_int(top);
	int2str(c_charge_value,&charge_value);

    for (int ind = 0; ind < in_para->size_population; ind++){
    	show_individual_msg(ind+1);
    	protein2pdbatoms(pdb_atoms, population_p[ind], top, population_p[ind]->z_matrix);
    	int2str(aux,&ind);
    	build_pdb_file_name(pdbfile_aux, aux, PREFIX_PDB_FILE_NAME_EA);
   	    save_pdb_file(in_para->path_local_execute, pdbfile_aux,&top->numatom,
    	    		pdb_atoms,NULL);
   	    compute_sas_multi(population_p[ind], in_para->path_local_execute,
   	    		in_para->path_gromacs_programs, pdbfile_aux, opt_fitness,
   	    		in_para->path_program_g_sas,in_para->path_program_read_g_sas,
   	    		in_para->computed_areas_g_sas_file,
   	    		&in_para->number_fitness,
   	    		in_para->path_program_clean_simulation);
   	    compute_energies_in_gromacs(population_p[ind],
   	    		in_para->path_local_execute, in_para->path_gromacs_programs,
   	    		pdbfile_aux,opt_fitness,c_charge_value,
   	    		in_para->path_program_compute_energy,
   	    		in_para->computed_energies_gromacs_file,
   	    		in_para->energy_file_xvg,in_para->path_program_read_energy,
   	    		in_para->computed_energy_value_file,
   	    		in_para->path_program_g_energy,
   	    		&in_para->number_fitness,
   	    		in_para->path_program_clean_simulation);
   	   compute_gyrate_mult(population_p[ind], in_para->path_local_execute,
   	    		in_para->path_gromacs_programs, pdbfile_aux, opt_fitness,
   	    		in_para->path_program_g_gyrate,
   	    		in_para->path_program_read_g_gyrate,
   	    		in_para->computed_radius_g_gyrate_file,
   	    		&in_para->number_fitness,
   	    		in_para->path_program_clean_simulation);
   	   compute_hbond_mult(population_p[ind], in_para->path_local_execute,
   	    		in_para->path_gromacs_programs, pdbfile_aux, opt_fitness,
   	    		in_para->path_program_g_hbond,
   	    		in_para->path_program_read_g_hbond,
   	    		&in_para->number_fitness,
   	    		in_para->path_program_clean_simulation);
   	compute_stride_mult(population_p[ind], in_para->path_local_execute,
	    		in_para->path_gromacs_programs, pdbfile_aux, opt_fitness,
	    		in_para->path_program_stride,
	    		&in_para->number_fitness,
	    		in_para->path_program_clean_simulation);



    }
    free(opt_fitness);
    desAllocate_pdbatom(pdb_atoms);
}


void build_random_individuals_multi(protein **p_new, const int *sizepop,
		const int *number_ind_select_reproduce){
	/* This function builds random individuals (Protein)
	 * how_many represents the number of individuals. Its value is obtained by
	 * difference between size population and number of individuals were
	 * selected to reproduction.
	 * last_index_individual_inserted represents the last index of individual
	 * which was inserted. In this moment, the new population is constituted
	 * by individuals which were selected to reproduction. Therefore, the last
	 * index value is *number_ind_select_reproduce -1 because the population
	 * index is started with 0.
	 */
	int how_many, last_index_individual_inserted;

	how_many = *sizepop - *number_ind_select_reproduce;
	last_index_individual_inserted = *number_ind_select_reproduce -1;

	build_random_individuals(p_new, &last_index_individual_inserted, &how_many);
}
