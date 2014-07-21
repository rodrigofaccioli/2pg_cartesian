#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "populationio.h"
#include "defines.h"
#include "futil.h"
#include "math_owner.h"
#include "pdbatom.h"
#include "topology.h"
#include "pdbio.h"
#include "pdbatom.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
void load_initial_population_file(protein_t *pop, const int *pop_size, const char *path, 
	const char *file_name, const primary_seq_t *primary_sequence){

	pdb_atom_t **atoms = NULL;
	const pdb_atom_t *atm_aux = NULL;
    int *num_atoms_by_model_PDB = NULL;
    char *path_pdb_file_name = NULL;


    //Loading PDB File of initial population
    path_pdb_file_name = path_join_file(path, 
        file_name);
    num_atoms_by_model_PDB = Malloc(int, *pop_size);
    get_num_atoms_by_model(num_atoms_by_model_PDB, path_pdb_file_name);

    atoms = allocate_Population_pdb(pop_size, num_atoms_by_model_PDB);
	load_pdb_model_file(atoms, NULL, path, file_name, num_atoms_by_model_PDB);

	//Allocation and Copy values to pop
	for (int i = 0; i < *pop_size; i++){
		pop[i].p_atoms = allocate_pdbatom(&num_atoms_by_model_PDB[i]);
		pop[i].p_topol = allocateTop_Global(&primary_sequence->num_res, 
			&num_atoms_by_model_PDB[i]);
		atm_aux = atoms[i];
		for (int a = 0; a < num_atoms_by_model_PDB[i]; a++){			
			copy_pdb_atom(&pop[i].p_atoms[a], &atm_aux[a]);
		}
		//Rename C-Terminal Oxygen atoms
		rename_oxygen_c_terminal(pop[i].p_atoms, &primary_sequence->num_res, 
			&num_atoms_by_model_PDB[i]);
	}
	//Building Topology of population
	build_topology_population(pop, pop_size);	

	free(num_atoms_by_model_PDB);
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
	
	for (int i  = 0; i < *num_model; i++){
		m = m + 1;
        if (m == 1){
            pdbfile = open_file(fname, fWRITE);
            writeHeader(pdbfile, 0.00, &pop[i].p_topol->numatom);
        }else{
            pdbfile = open_file(fname, fAPPEND);
        }		
		writeModel(pdbfile, &m);
		writeATOM(pdbfile, pop[i].p_atoms, &pop[i].p_topol->numatom);
		writeEndModel(pdbfile);
		fclose(pdbfile);
	}	
	free(fname);	
}
