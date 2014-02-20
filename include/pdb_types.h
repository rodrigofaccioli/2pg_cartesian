#ifndef OLD_PDB_TYES_H
#define OLD_PDB_TYES_H

#define NUM_ATOM_NAME 5
#define NUM_RES_NAME 4

#include "enums.h"
#include "vector_types.h"

typedef struct spdbatom{
	type_atoms_t atomid;
    char atmname[NUM_ATOM_NAME];
    char resname[NUM_RES_NAME];
    int resnum;
    own_vector_t coord;
    int atmnumber;
 }pdb_atom_t;

 typedef struct //pdbseqres
  {
    char chain[1];
    char resname[3];
  }pdb_seqres_t;

#endif