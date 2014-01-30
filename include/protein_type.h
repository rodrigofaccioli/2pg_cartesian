#ifndef OLD_PROTEIN_TYPE_H
#define OLD_PROTEIN_TYPE_H

#include "topology_types.h"
#include "pdb_types.h"

/**
 * protein_t represents an atomistic conformation of protein
 * ID is an identification of protein
 * p_atoms pointer to atoms of protein
 * p_topol pointer to topology of protein 
*/
typedef struct sprotein {
   int ID;
   pdb_atom_t *p_atoms;
   top_global_t *p_topol;
 }protein_t;


#endif
