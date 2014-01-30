#include "topologylib.h"
#include "aminoacids_io.h"

/** _load_amino_seq loads a Fasta File

*/

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
    char line[MAX_LINE_FASTA+1]; //represents Fasta lines
    char seq_line[MAX_LEN_PROTEIN];//represents the primary sequence of protein
    char c;
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
	remove_character(line, '\n');
	strcpy(seq_line,line);
	while ( (fgets(line,MAX_LINE_FASTA,arq) != NULL) &&
			(read_fasta_file == btrue) ){
		remove_character(line,'\n');
		trim(line);
		if (strncmp(line,">",1) == 0){ //Force to work with one sequence only
			read_fasta_file = bfalse;
		}else{
			strcat(seq_line,line);
		}
	}
	fclose(arq);
	remove_character(seq_line, '\n');
	n_residues = strlen(seq_line);
	seq_prim = allocate_primary_seq(&n_residues){
	
	for (i = 0; i < *n_residues; i++){
		c = seq_line[i];
		strcpy(seq_prim->seq_res[i].id_1, c);
		seq_prim->seq_res[i].id = get_amino_id(c);

        seq_index = seq_index + 1;
	}
	return seq_prim;

}

