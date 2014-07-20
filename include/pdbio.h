#ifndef OLD_PDBIO_H
#define OLD_PDBIO_H

#include <stdio.h>
#include "pdbatom.h"

#define MAX_LINE_PDB 100

_2PG_CARTESIAN_EXPORT
void save_pdb_file(const char *path, const char *file_name, const int *numatom,
		const pdb_atom_t *atoms, const pdb_seqres_t *seqres);
_2PG_CARTESIAN_EXPORT
void show_coordinates(const pdb_atom_t *atoms, const int *numatom);
_2PG_CARTESIAN_EXPORT
void show_coordinate(const pdb_atom_t *atom, const int *atom_index);
_2PG_CARTESIAN_EXPORT
void load_pdb_file(pdb_atom_t *atoms, pdb_seqres_t *seqres,
		const char *path, const char *file_name, const int *numatom);
_2PG_CARTESIAN_EXPORT
void load_pdb_file_without_num_atom(pdb_atom_t *atoms, pdb_seqres_t *seqres,
		const char *path_file_name);
_2PG_CARTESIAN_EXPORT
void load_pdb_model_file(pdb_atom_t **atoms, pdb_seqres_t *seqres,
		const char *path, const char *file_name, const int *num_atoms_by_model);
_2PG_CARTESIAN_EXPORT
int get_num_atom(const char *path_PDB_file_name);
_2PG_CARTESIAN_EXPORT
void get_num_atoms_by_model(int *num_atoms_by_model, const char *path_PDB_file_name);
_2PG_CARTESIAN_EXPORT
void save_model_pdb_file(const char *path, const char *file_name, const int *num_model, 
	const int *numatom, pdb_atom_t **atoms_model, const pdb_seqres_t *seqres );
_2PG_CARTESIAN_EXPORT
void load_pdb_atoms(char line[], pdb_atom_t *atoms, const int *l);

_2PG_CARTESIAN_EXPORT
void writeHeader(FILE *pdbfile, float dif, const int *npdb );
_2PG_CARTESIAN_EXPORT
void writeATOM(FILE *pdbfile, const pdb_atom_t *atoms, const int *npdb );
_2PG_CARTESIAN_EXPORT
void writeModel(FILE *pdbfile, const int *model);
_2PG_CARTESIAN_EXPORT
void writeEndModel(FILE *pdbfile);
#endif


