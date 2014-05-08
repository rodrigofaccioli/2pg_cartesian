#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "defines.h"
#include "enums.h"
#include "protein.h"
#include "topology.h"
#include "pdbio.h"
#include "pdbatom.h"
#include "messages.h"
#include "algorithms.h"
#include "string_owner.h"
#include "futil.h"
#include "random_number_gsl.h"
#include "aminoacids.h"
#include "aminoacids_io.h"
#include "populationio.h"
#include "topology.h"
#include "topologyio.h"
#include "rotation.h"
#include "math_owner.h"
#include "solution.h"
#include "gromacs_objectives.h"
#include "algorithms.h"
#include "randomlib.h"
#include "objective.h"
#include "solutionio.h"

void save_energies(const solution_t *solution_curr, 
    const solution_t *solution_new, const input_parameters_t *in_para, 
    const int *num_sol){    
    const protein_t* p_curr;
    FILE * energy_file; 
    char *file_name;
    file_name = Malloc(char, MAX_FILE_NAME);
    strcpy(file_name, "monte_carlo_energies.fit");
    char *fname = path_join_file(in_para->path_local_execute,
        file_name);    
    p_curr = (protein_t*) solution_curr->representation;
    if (*num_sol == 1){
        energy_file = open_file(fname, fWRITE);
        fprintf (energy_file,"#Index New_Solution \t Current_Solution \n");
    }else{
        energy_file = open_file(fname, fAPPEND);
    }    
    fprintf (energy_file,"%i \t %f \t %f \n", *num_sol, solution_new[0].obj_values[0], solution_curr[0].obj_values[0]);
    fclose(energy_file);    
    free(fname);
    free(file_name);    
}

void save_solution(const solution_t *solution_curr, 
    const input_parameters_t *in_para, const int *model){
    const protein_t* p_curr;
    FILE * pdbfile; 
    char *file_name;
    file_name = Malloc(char, MAX_FILE_NAME);
    strcpy(file_name, "monte_carlo_solutions.pdb");
    char *fname = path_join_file(in_para->path_local_execute,
        file_name);    
    p_curr = (protein_t*) solution_curr->representation;
    if (*model == 1){
        pdbfile = open_file(fname, fWRITE);
        writeHeader(pdbfile, 0.00, &p_curr->p_topol->numatom);
    }else{
        pdbfile = open_file(fname, fAPPEND);
    }
    writeModel(pdbfile, model);
    writeATOM(pdbfile, p_curr->p_atoms, &p_curr->p_topol->numatom);
    writeEndModel(pdbfile);
    fclose(pdbfile);

    free(fname);
    free(file_name);    
}

/** Accepts updating the current solution with new solution
* solution_curr is the current solution
* solution_new is the new solution
*/ 
void accept(solution_t *solution_curr, const solution_t *solution_new){
	protein_t *p_curr;
	const protein_t *p_new;
	p_curr = (protein_t*) solution_curr[0].representation;
	p_new = (protein_t*) solution_new[0].representation;
	copy_protein(p_curr, p_new);
	solution_curr->obj_values[0] = solution_new->obj_values[0];
}

/** Updates new solution with current solution
* solution_new is the new solution
* solution_curr is the current solution
*/ 
void reject(solution_t *solution_new, const solution_t *solution_curr){
	const protein_t *p_curr;
	protein_t *p_new;
	p_curr = (protein_t*) solution_curr[0].representation;
	p_new = (protein_t*) solution_new[0].representation;
	copy_protein(p_new, p_curr);
	solution_new->obj_values[0] = solution_curr->obj_values[0];
}


static int compare_solution(const void *x, const void *y){
    double fx, fy;
    fx = ((solution_t *)x)->obj_values[0];
    fy = ((solution_t *)y)->obj_values[0];
    if (fx > fy){
        return 1;
    }else if ( (fx -fy) == 0.0001){
        return 0;
    }else{
        return -1;
    }
}

/** Builds a new solution and computes the objective
* ind_new is an protein structure
* in_para is the input parameter
*/
void mc_build_solution(protein_t *ind_new, const input_parameters_t *in_para){
    float angle_radians, angle_degree, rate;
    int what_rotation;
    int max_kind_of_rotation = 4;
    int num_residue_choose;
    
    //Choose a residue. It must be started 1 untill number of residue
    num_residue_choose = get_choose_residue(&ind_new->p_topol->numres);
    //Obtaing kind of rotation
    what_rotation = _get_int_random_number(&max_kind_of_rotation);            
    //Appling the rotation 
    if (what_rotation == 0){
        //Obtaing a random degree angule 
        angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_psi, 
        &in_para->max_angle_mutation_psi);
        //Cast to radians
        angle_radians = degree2radians(&angle_degree);
        rotation_psi(ind_new, &num_residue_choose, &angle_radians);
    }else if (what_rotation == 1){
        //Obtaing a random degree angule 
        angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_phi, 
        &in_para->max_angle_mutation_phi);
        //Cast to radians
        angle_radians = degree2radians(&angle_degree);
        rotation_phi(ind_new, &num_residue_choose, &angle_radians);            
    }else if (what_rotation == 2){
        //Obtaing a random degree angule 
        angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_omega, 
        &in_para->max_angle_mutation_omega);
        //Cast to radians
        angle_radians = degree2radians(&angle_degree);
        rotation_omega(ind_new, &num_residue_choose, &angle_radians);            
    }else if (what_rotation == 3){
        int chi = 0;
        int max_chi;
        int num_chi_choose;
        char *res_name;
        res_name = Malloc(char, 4);
        get_res_name_from_res_num(res_name, &num_residue_choose, 
            ind_new->p_atoms, &ind_new->p_topol->numatom);
        max_chi = get_number_chi(res_name);
        if (max_chi > 0){
            //Choose a chi of residue. It must be started 1 untill number of residue
            if (max_chi == 1){
                chi = 1;
            }else{
                while (chi == 0){                
                    num_chi_choose = max_chi + 1;
                    chi = _get_int_random_number(&num_chi_choose);
                }                    
            }                
            //Obtaing a random degree angule 
            angle_degree = _get_float_random_interval(&in_para->min_angle_mutation_side_chain, 
            &in_para->max_angle_mutation_side_chain);
            //Cast to radians
            angle_radians = degree2radians(&angle_degree);
            //printf("what_rotation %i num_residue_choose %i chi %i \n", 3, num_residue_choose, chi);
            rotation_chi(ind_new, &num_residue_choose, &chi, &angle_radians);
        }
        free(res_name);
    }    
}

int mc_metropolis(const input_parameters_t *in_para){

    char *prefix;
    int num_solution; // stores the number of solution
    primary_seq_t *primary_sequence; // Primary Sequence of Protein
    
    protein_t *prot_new; // new protein structure. It is build by mc_build_solution function
    protein_t *prot_curr; // current protein structure

    solution_t *solution_new; // new solution
    solution_t *solution_curr; // current solution
    float R;
    float prob, rr;
    int T;
    int model;

    //Starting values
    num_solution = 1;
    rr = 0;
	R = 8.3144621E-3; // kJ/Kmol
	T = 309;	// K	
	prob = 1.000;
    model = 0;

    //Loading Fasta file
    primary_sequence = _load_amino_seq(in_para->seq_protein_file_name);

    //Allocating PDB ATOMS
    prot_new = allocateProtein(&num_solution);
    prot_curr = allocateProtein(&num_solution);
    
    //Loading initial population
    load_initial_population_file(prot_curr, &num_solution, 
        in_para->path_local_execute,in_para->initial_pop_file_name,
        primary_sequence);    

    //Setting protein previous that is the same of prot_curr
    copy_protein(&prot_new[0], &prot_curr[0]);

    //Saving topology of population 
    prefix = Malloc(char,10);
    strcpy(prefix, "prot");
    save_topology_population(prot_curr, &num_solution, 
    in_para->path_local_execute, prefix);    
    free(prefix);

/**************** STARTING Monte Carlo Metropolis Algorithm *************************/
    display_msg("Starting Monte Carlo Metropolis Algorithm \n");
    initialize_algorithm_execution(primary_sequence, in_para);
    init_gromacs_execution();    
    //Setting solutions
    solution_curr = allocate_solution(&num_solution, &in_para->number_fitness);
    solution_new = allocate_solution(&num_solution, &in_para->number_fitness);    
    //Setting reference of proteins to solution 
    set_proteins2solutions(solution_curr, prot_curr, &num_solution);
    set_proteins2solutions(solution_new, prot_new, &num_solution);
	//Computing solutions with GROMACS
    get_gromacs_objectives(solution_curr, in_para);
    for (int s = 1; s <= in_para->MonteCarloSteps; s++){
        //Saving PDB of current solution
        if (s % in_para->freq_mc == 0){
            model  = model + 1;
            save_solution(&solution_curr[0], in_para, &model);
        }
        //Building new Solution
    	mc_build_solution(&prot_new[0], in_para);
        //Calculates energy of new Solution
    	get_gromacs_objectives(&solution_new[0], in_para);
        //Saving energies of New and Current solution
        save_energies(solution_curr, solution_new, in_para, &s);
		// Checking the new solution acceptance		
		if( solution_new[0].obj_values[0] > solution_curr[0].obj_values[0] ){// If the energy of the new structure is higher than of the previous...
			prob=exp( (-1)*((solution_new[0].obj_values[0] - solution_curr[0].obj_values[0])/(R*T) ) );// Metropolis criterion
			rr = _get_double_random_number();
            printf("prob %f rr %f \n", prob, rr);
			if( rr <= prob ){
				accept(&solution_curr[0],&solution_new[0]);
			}else{
				reject(&solution_new[0],&solution_curr[0]);
			}
		}else{ 
			// new is better than current
			accept(&solution_curr[0], &solution_new[0]);
		}
    }    
    //save_final_solution(&solution_curr[0], in_para);
    finish_gromacs_execution();

/**************** FINISHED Mono-objetive Evolutionary Algorithm *************************/     

    desallocate_solution(solution_curr, &num_solution);
    desallocate_solution(solution_new, &num_solution);
    desallocateProtein(prot_curr, &num_solution);
    desallocateProtein(prot_new, &num_solution);
    desallocate_primary_seq(primary_sequence);
    _finish_random_gsl();

	return 0;
}
