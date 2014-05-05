#ifndef OLD_PDBIO_H
#define OLD_PDBIO_H

#include <stdio.h>
#include "pdbatom.h"

#define MAX_LINE_PDB 100

void save_pdb_file(const char *path, const char *file_name, const int *numatom,
		const pdb_atom_t *atoms, const pdb_seqres_t *seqres);
void show_coordinates(const pdb_atom_t *atoms, const int *numatom);
void show_coordinate(const pdb_atom_t *atom, const int *atom_index);
void load_pdb_file(pdb_atom_t *atoms, pdb_seqres_t *seqres,
		const char *path, const char *file_name, const int *numatom);
void load_pdb_file_without_num_atom(pdb_atom_t *atoms, pdb_seqres_t *seqres,
		const char *path_file_name);
void load_pdb_model_file(pdb_atom_t **atoms, pdb_seqres_t *seqres,
		const char *path, const char *file_name, const int *num_atoms_by_model);
int get_num_atom(const char *path_PDB_file_name);
void get_num_atoms_by_model(int *num_atoms_by_model, const char *path_PDB_file_name);
void save_model_pdb_file(const char *path, const char *file_name, const int *num_model, 
	const int *numatom, pdb_atom_t **atoms_model, const pdb_seqres_t *seqres );
void load_pdb_atoms(char line[], pdb_atom_t *atoms, const int *l);

void writeHeader(FILE *pdbfile, float dif, const int *npdb );
void writeATOM(FILE *pdbfile, const pdb_atom_t *atoms, const int *npdb );
void writeModel(FILE *pdbfile, const int *model);
void writeEndModel(FILE *pdbfile);
#endif


