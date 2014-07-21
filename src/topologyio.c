#include<stdio.h>
#include<stdlib.h>
#include<string.h>


#include "topologyio.h"
#include "defines.h"
#include "futil.h"
#include "enums.h"
#include "messages.h"
#include "string_owner.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

#define MAX_LINE_FASTA 81

static void build_topol_filename(char *filename, const char *prefix, const int *v){

	char *aux;
	aux = Malloc(char, 10);

	int2str(aux,v);
	strcpy(filename, prefix);
	strcat(filename,"_");
	strcat(filename,aux);
	strcat(filename,".top");
	free(aux);
}

static void write_header_top_residue_atom_info_psi(FILE * topolfile){
	char header [] = "; Atoms of the PSI rotation";
	fprintf(topolfile,"%s\n", header);
}

static void write_header_top_residue_atom_info_phi(FILE * topolfile){
	char header [] = "; Atoms of the PHI rotation";
	fprintf(topolfile,"%s\n", header);
}

static void write_header_top_residue_atom_info_omega(FILE * topolfile){
	char header [] = "; Atoms of the OMEGA rotation";
	fprintf(topolfile,"%s\n", header);
}

static void write_header_top_residue_atom_info_side_chains(FILE * topolfile){
	char header [] = "; Atoms of the Side Chains rotation";
	fprintf(topolfile,"%s\n", header);
}

static void write_top_residue_atom_info_withoutmoved(FILE *topolfile, const int *res,
	const top_residue_atom_info_t *info){
	char *line, *aux;

	line = Malloc(char, MAX_LINE_FILE);
	aux = Malloc(char, 10);
	fprintf(topolfile,"\tResidue: %i\n",*res+1);
	//Printing fixed atoms
	strcpy(line, "\tFixed Atoms: ");
	for (int n = 0; n < info[*res].num_fixed; n++){
		int2str(aux, &info[*res].fixed_atoms[n]);
		strcat(line, aux);
		strcat(line, " ");
	}
	fprintf(topolfile,"%s\n\n",line);
	
	free(aux);
	free(line);
}

static void write_top_residue_atom_info(FILE *topolfile, const int *res,
	const top_residue_atom_info_t *info){
	char *line, *aux;

	line = Malloc(char, MAX_LINE_FILE);
	aux = Malloc(char, 10);
	fprintf(topolfile,"\tResidue: %i\n",*res+1);
	//Printing fixed atoms
	strcpy(line, "\tFixed Atoms: ");
	for (int n = 0; n < info[*res].num_fixed; n++){
		int2str(aux, &info[*res].fixed_atoms[n]);
		strcat(line, aux);
		strcat(line, " ");
	}
	fprintf(topolfile,"%s\n",line);
	//Printing moved atoms
	strcpy(line, "\tMoved Atoms: ");
	for (int n = 0; n < info[*res].num_moved; n++){
		int2str(aux, &info[*res].moved_atoms[n]);
		strcat(line, aux);
		strcat(line, " ");
	}
	fprintf(topolfile,"%s\n",line);
	fprintf(topolfile,"\n");

	free(aux);
	free(line);
}

static void write_top_residue_atom_info_side_chain(FILE *topolfile, const int *res,
	const top_residue_side_chains_t *side_chains){
	char *line, *aux;
		 
	line = Malloc(char, MAX_LINE_FILE);
	aux = Malloc(char, 10);
	fprintf(topolfile,"\tResidue: %i\n",*res+1);
	if (side_chains[*res].num_chi == 0){
		fprintf(topolfile,"\tNO SIDE CHAIN\n\n");
	}else{
		for (int chi = 1; chi <= side_chains[*res].num_chi;chi++){
			strcpy(line, "\tchi ");
			int2str(aux, &chi);
			strcat(line, aux);
			strcat(line, ":\n");
			fprintf(topolfile, line);
			//Printing fixed atoms
			strcpy(line, "\tFixed Atoms: ");
			for (int n = 0; n < side_chains[*res].atoms_chi[chi-1].num_fixed; n++){
				int2str(aux, &side_chains[*res].atoms_chi[chi-1].fixed_atoms[n]);
				strcat(line, aux);
				strcat(line, " ");
			}
			fprintf(topolfile,"%s\n",line);
			//Printing moved atoms
			strcpy(line, "\tMoved Atoms: ");
			for (int n = 0; n < side_chains[*res].atoms_chi[chi-1].num_moved; n++){
				int2str(aux, &side_chains[*res].atoms_chi[chi-1].moved_atoms[n]);
				strcat(line, aux);
				strcat(line, " ");
			}
			fprintf(topolfile,"%s\n",line);			
		}
		fprintf(topolfile,"\n");
	}	
	free(aux);
	free(line);
}


static void write_top_residue_atom_info_psi(FILE *topolfile, const top_residue_atom_info_t *psi, 
	const int *num_res){
	write_header_top_residue_atom_info_psi(topolfile); 
	for (int r = 0; r < *num_res-1; r++){
		write_top_residue_atom_info(topolfile, &r, psi);
	}	
}

static void write_top_residue_atom_info_phi(FILE *topolfile, const top_residue_atom_info_t *phi, 
	const int *num_res){
	write_header_top_residue_atom_info_phi(topolfile); 
	for (int r = 1; r < *num_res; r++){
		write_top_residue_atom_info(topolfile, &r, phi);
	}	
}

static void write_top_residue_atom_info_omega(FILE *topolfile, const top_residue_atom_info_t *omega, 
	const int *num_res){
	write_header_top_residue_atom_info_omega(topolfile); 
	for (int r = 0; r < *num_res-1; r++){
		write_top_residue_atom_info_withoutmoved(topolfile, &r, omega);
	}	
}

static void write_top_residue_atom_info_side_chains(FILE *topolfile, 
	const top_residue_side_chains_t *side_chains, const int *num_res){
	write_header_top_residue_atom_info_side_chains(topolfile);
	for (int r = 0; r < *num_res; r++){
		write_top_residue_atom_info_side_chain(topolfile, &r, side_chains);
	}	
}

static void write_top_residue_range_atoms(FILE *topolfile, const top_residue_range_atoms_t *range, 
	const int *num_res){
	char topic [] = "; Residues and Its Range Atoms";
	fprintf(topolfile,"%s\n", topic);
	char header [] = "; Residue First_Atom Last_Atom ";
	fprintf(topolfile,"%s\n", header);
	
	for (int r = 0; r < *num_res; r++){
		fprintf(topolfile,"%i %i %i\n", r+1, range[r].first_atom, range[r].last_atom);	
	}
	fprintf(topolfile,"\n \n");
}

void save_topology_protein(const top_global_t *top, const char *path, 
	const char *filename){
	FILE *topolfile=NULL;	
	char *fname = path_join_file(path,filename);
	topolfile = open_file(fname, fWRITE);
	write_top_residue_range_atoms(topolfile, top->range_atoms, &top->numres);
	write_top_residue_atom_info_psi(topolfile, top->psi, &top->numres);
	write_top_residue_atom_info_phi(topolfile, top->phi, &top->numres);
	write_top_residue_atom_info_omega(topolfile, top->omega, &top->numres);
	write_top_residue_atom_info_side_chains(topolfile, top->side_chains, &top->numres);

	fclose(topolfile);
	free(fname);
}

/** Save topology of population
*/
_2PG_CARTESIAN_EXPORT
void save_topology_population(const protein_t *pop, const int *popsize, 
	const char *path, const char *prefix){
	char *filename;
	filename = Malloc(char, MAX_FILE_NAME);
	for (int i = 0; i < *popsize; i++){
		build_topol_filename(filename, prefix, &i);
		save_topology_protein(pop[i].p_topol, path, filename);
	}
	free(filename);

}

static boolean_t _check_pdb_fasta_file(char *line){
	/*This function check if line contains pdbid which means that
	 * the file is pdb fasta file
	 * Example line: >1VII:A|PDBID|CHAIN|SEQUENCE
	 */
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

void _create_fasta_pdb(const char *prot_name, const char *chain_name,
		const char *prot_seq, const char *file_name_protein){
	/* Create a pdb fasta file.
	 * Example of pdb fasta file:
	 * >1VII:A|PDBID|CHAIN|SEQUENCE
     * MLSDEDFKAVFGMTRSAFANLPLWKQQNLKKEKGLF
	*/
	FILE *arq;
	int seq_index = -1;
	char *prot_seq_aux;
	int len = 0;
	char header[MAX_LINE_FASTA];

	arq = open_file(file_name_protein, fWRITE);
	//Writing header
	strcpy(header,">");
	strcat(header,prot_name);
	strcat(header,":");
	strcat(header,chain_name);
	strcat(header,"|PDBID|CHAIN|SEQUENCE");
	fprintf (arq, "%s\n", header);
	len = strlen(prot_seq);
	prot_seq_aux = Malloc(char, MAX_LINE_FASTA+1);

	strncpy(prot_seq_aux,prot_seq,MAX_LINE_FASTA);
	fprintf (arq, "%s\n",prot_seq_aux);
    free(prot_seq_aux);
	fclose(arq);

}
