#include<stdio.h>
#include<string.h>

#include "defines.h"
#include "enums.h"
#include "topologylib.h"
#include "aminoacids_io.h"
#include "protein.h"
#include "futil.h"
#include "string_owner.h"
#include "messages.h"
#include "aminoacids.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

#define MAX_LINE_FASTA 81
#define MAX_LEN_PROTEIN 1000


/** _check_pdb_fasta_file checks if line contains a pdbid which means that
* this file is pdb fasta.
* Example line: >1VII:A|PDBID|CHAIN|SEQUENCE
*/
static boolean_t _check_pdb_fasta_file(char *line){
	char *pch;
	pch = strtok (line,"|");
	while (pch != NULL){
		if (strcmp(pch,"PDBID") == 0){
			return btrue;
		}
		pch = strtok (NULL, "|");
	}
	return bfalse;
}


/** _load_amino_seq loads a Fasta File

*/
_2PG_CARTESIAN_EXPORT
primary_seq_t *_load_amino_seq(const char *file_name_protein){	
	/* load amino_t based on pdb fasta file.
	 * Example of pdb fasta file:
	 * >1VII:A|PDBID|CHAIN|SEQUENCE
     * MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF
	*/
	FILE *arq;
	primary_seq_t * seq_prim;
	int n_residues;
	int fscanfError;	
    int seq_index;    
    char line[MAX_LINE_FASTA+1] = "\0"; //represents Fasta lines
    char seq_line[MAX_LEN_PROTEIN] = "\0";//represents the primary sequence of protein    
    char aux_line;    
    int i;
    boolean_t first_line = btrue;
    boolean_t read_fasta_file = btrue;

	arq = open_file(file_name_protein, fREAD);
	fgets(line,MAX_LINE_FASTA,arq);
	if (_check_pdb_fasta_file(line) == bfalse){
		char msg[300];
	    sprintf(msg,"Fasta file %s is not a pdb fasta file \n",file_name_protein);
		fatal_error(msg);
	}
	//After checking Fasta file, initialize variables
	seq_index = 0;
	/*Getting the primary sequence. It is represented by seq_line variable.
	 * So it contains residues from one chain of protein only. Because of this
	 * is checked if line starts with >*/
	fgets(line,MAX_LINE_FASTA,arq); //Here line contains residues of protein. Maybe it is first part like 1E8A
	trim(line);
	remove_character_enter(line);
	strcpy(seq_line,line);
	while ( (fgets(line,MAX_LINE_FASTA,arq) != NULL) &&
			(read_fasta_file == btrue) ){
		remove_character_enter(line);
		trim(line);
		if (strncmp(line,">",1) == 0){ //Force to work with one sequence only
			read_fasta_file = bfalse;
		}else{
			strcat(seq_line,line);
		}
	}
	fclose(arq);
	remove_character_enter(seq_line);
	n_residues = strlen(seq_line);
	seq_prim = allocate_primary_seq(&n_residues);

	check_terminal_charge(seq_line, &n_residues);
	for (int i = 0; i < n_residues; i++){
		aux_line = seq_line[i];
		strcpy(seq_prim->seq_res[i].id_1, &aux_line);
		seq_prim->seq_res[i].id = _get_amino_id_1(aux_line);
		set_amino_id_3(seq_prim->seq_res[i].id_3, &seq_prim->seq_res[i].id);
        seq_index = seq_index + 1;
	}
	return seq_prim;

}
