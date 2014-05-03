#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "populationio.h"
#include "futil.h"
#include "math_owner.h"
#include "pdbatom.h"
#include "topology.h"
#include "pdbio.h"
#include "pdbatom.h"


void load_initial_population_file(protein_t *pop, const int *pop_size, const char *path, 
	const char *file_name, const primary_seq_t *primary_sequence){

	pdb_atom_t **atoms;
	const pdb_atom_t *atm_aux;
    int num_atoms_PDB;
    char *path_pdb_file_name;


    //Loading PDB File of initial population
    path_pdb_file_name = path_join_file(path, 
        file_name);
    num_atoms_PDB  = get_num_atom(path_pdb_file_name);
    atoms = allocate_Population_pdb(pop_size, &num_atoms_PDB);
	load_pdb_model_file(atoms,NULL, path, file_name, &num_atoms_PDB);

	//Allocation and Copy values to pop
	for (int i = 0; i < *pop_size; i++){
		pop[i].p_atoms = allocate_pdbatom(&num_atoms_PDB);
		pop[i].p_topol = allocateTop_Global(&primary_sequence->num_res, 
			&num_atoms_PDB);
		atm_aux = atoms[i];
		for (int a = 0; a < num_atoms_PDB; a++){			
			copy_pdb_atom(&pop[i].p_atoms[a], &atm_aux[a]);
		}
		//Rename C-Terminal Oxygen atoms
		rename_oxygen_c_terminal(pop[i].p_atoms, &primary_sequence->num_res, 
			&num_atoms_PDB);
	}

	//Building Topology of population
	build_topology_population(pop, pop_size);

	free(path_pdb_file_name);
	desAllocate_Population_pdb(atoms, pop_size);

}

/** Save a set of models in PDB format
* pop is population
* path is the path where file_name will  be saved
* file_name is the name of file
* num_model the number of models. It can be be the number of individuals
*/
void save_population_file(const protein_t *pop, const char *path, const char *file_name, 
	const int *num_model ){
	FILE *pdbfile=NULL;
	int m = 0;
	char *fname = path_join_file(path,file_name);
	pdbfile = open_file(fname, fWRITE);
	writeHeader(pdbfile, 10.5, &pop[0].p_topol->numatom);
	for (int i  =0; i < *num_model; i++){
		m = m + 1;
		writeModel(pdbfile, &m);
		writeATOM(pdbfile, pop[i].p_atoms, &pop[i].p_topol->numatom);
		writeEndModel(pdbfile);
	}	
	free(fname);
	fclose(pdbfile);
}