#ifndef OLD_AMINOACIDS_TYPES_H
#define OLD_AMINOACIDS_TYPES_H

#include "enums.h"

/** amino_t represents one aminoacid (residue)
*/

#ifdef __cplusplus
extern "C"
{
#endif

typedef struct samino {
   type_aminos_t id;   
   char *id_1; // code for 1 Letter
   char *id_3; // code for 3 Letters
 }amino_t;


/** primary_seq_t represents the primary sequence of protein
*/
typedef struct sprimary_seq {
   int num_res;  // number of residues
   amino_t * seq_res; // sequence of residues
 }primary_seq_t;

#ifdef __cplusplus
}
#endif

 #endif


 