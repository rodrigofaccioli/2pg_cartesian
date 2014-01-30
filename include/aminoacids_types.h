/** amino_t represents one aminoacid (residue)
*/

typedef struct samino {
   type_aminos_t id;   
   char id_1[2]; // code for 1 Letter
   char id_3[4]; // code for 3 Letters
 }amino_t;


/** primary_seq_t represents the primary sequence of protein
*/
typedef struct sprimary_seq {
   int num_res;  // number of residues
   amino_t * seq_res; // sequence of residues
 }primary_seq_t;


 