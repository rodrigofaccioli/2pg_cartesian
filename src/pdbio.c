#include<stdio.h>
#include<string.h>

#include "defines.h"
#include "enums.h"
#include "protein.h"
#include "pdbio.h"
#include "pdbatom.h"
#include "futil.h"
#include "messages.h"
#include "string_owner.h"

void load_pdb_atoms(char line[], pdb_atom_t *atoms, const int *l){
	char *section;
	int index_atom;
	char *atmname;
	char *resname_aux;
	char *chain_name;
	int atm_number;
	char *atm_number_c;
	int resnum;
	char *resnum_c;
	char *coord_c;
	float x,y,z;
	//Allocating
	section = Malloc(char,7);
	atmname = Malloc(char,6);
	resname_aux = Malloc(char,4);
	chain_name = Malloc(char,2);
	atm_number_c = Malloc(char,6);
	resnum_c = Malloc(char,10);//Malloc(char,5);
	coord_c = Malloc(char,9);

	//Section
	substring(section,line,0,5);
	//Atom Number
	substring(atm_number_c,line,6,11);
	trim(atm_number_c);
	atm_number = str2int(atm_number_c);
	//Atom Name
	substring(atmname,line,11,16);
	trim(atmname);
	//ChainID
	//substring(chain_name,line,15,16);
	strcpy(chain_name," ");
	//Residue Name
	substring(resname_aux,line,17,20);
	trim(resname_aux);
	//Residue Number
	substring(resnum_c,line,22,28);
	trim(resnum_c);
	resnum = str2int(resnum_c);
	//Coordinates
	substring(coord_c,line,30,38);
	trim(coord_c);
	x = str2float(coord_c);
	substring(coord_c,line,39,46);
	trim(coord_c);
	y = str2float(coord_c);
	substring(coord_c,line,47,55);
	trim(coord_c);
	z = str2float(coord_c);

	set_pdb_atom_coordinates(&atoms[*l], atmname, resname_aux, chain_name, &resnum,
				&x, &y, &z,&atm_number);
	free(section);
	free(atmname);
	free(resname_aux);
	free(chain_name);
	free(atm_number_c);
	free(resnum_c);
	free(coord_c);

}

static void writeHeader(FILE *pdbfile, float dif, const int *npdb ){
	fprintf(pdbfile,"HEADER    PROTEIN %lf %d \n",dif,*npdb);
	fprintf(pdbfile,"COMPND\n");
	fprintf(pdbfile,"SOURCE\n");
}

static void writeModel(FILE *pdbfile, const int *model){
	fprintf(pdbfile,"MODEL        %i\n",*model);
}

static void writeEndModel(FILE *pdbfile){
	fprintf(pdbfile,"TER\n");
	fprintf(pdbfile,"ENDMDL\n");
}

static void writeATOM(FILE *pdbfile, const pdb_atom_t *atoms, const int *npdb ){
	char chnname[2];
	char resname_aux[4];
	strcpy(chnname, "A");
	for (int i=0; i<*npdb; i++)	{
             if ( strcmp(atoms[i].resname, "CYX") == 0)
		          strcpy (resname_aux, "CYS");
             else
            	 strcpy (resname_aux, atoms[i].resname);
             fprintf(pdbfile,"%6.6s%5d %4.4s %3.3s %s%4d    %8.3f%8.3f%8.3f \n","ATOM  ",i+1,atoms[i].atmname,resname_aux,
            		             chnname,atoms[i].resnum,atoms[i].coord.x,atoms[i].coord.y,atoms[i].coord.z);
	}
}


void save_pdb_file(const char *path, const char *file_name, const int *numatom,
		const pdb_atom_t *atoms, const pdb_seqres_t *seqres){
	FILE *pdbfile=NULL;
	char *fname = path_join_file(path,file_name);
	pdbfile = open_file(fname, fWRITE);
	writeHeader(pdbfile, 10.5, numatom);
	writeATOM(pdbfile, atoms, numatom);
	free(fname);
	fclose(pdbfile);
}

void save_model_pdb_file(const char *path, const char *file_name, const int *num_model, 
	const int *numatom, pdb_atom_t **atoms_model, const pdb_seqres_t *seqres ){
	/* This function save a set of models in PDB format.
	   atoms_model can be a population of pdb_atom_t
	*/
	FILE *pdbfile=NULL;
	int m = 0;
	char *fname = path_join_file(path,file_name);
	pdbfile = open_file(fname, fWRITE);
	writeHeader(pdbfile, 10.5, numatom);
	for (int i  =0; i < *num_model; i++){
		m = m + 1;
		writeModel(pdbfile, &m);
		writeATOM(pdbfile, atoms_model[i], numatom);
		writeEndModel(pdbfile);
	}	
	free(fname);
	fclose(pdbfile);
}

void show_coordinates(const pdb_atom_t *atoms, const int *numatom){
	/*Show coordinates of atoms of protein*/
	display_msg("Showing coordinates of atoms \n");
	display_msg("atom x y z \n");
	for (int i = 0; i < *numatom; i++){
		show_coordinate(&atoms[i],&i);
	}
}

void show_coordinate(const pdb_atom_t *atom, const int *atom_index){
	/*Show coordinate of an atom of protein
	 * If you want to see all coordinates, see show_coordinates function
	 */
	char msg[50];
	sprintf(msg,"%i %f %f %f \n", *atom_index, atom->coord.x, atom->coord.y,
			atom->coord.z);
	display_msg(msg);
}

void load_pdb_file(pdb_atom_t *atoms, pdb_seqres_t *seqres,
		const char *path, const char *file_name, const int *numatom){
	FILE *pdbfile=NULL;
	char pdb_line [MAX_LINE_PDB];
	char *aux = NULL;
	int line_atm = -1;
        char *fname = path_join_file(path,file_name);
	pdbfile = open_file(fname, fREAD);
	while ( fgets(pdb_line,MAX_LINE_PDB,pdbfile) != NULL){
		if (strncmp(pdb_line,"ATOM",4) == 0){
			line_atm = line_atm + 1;
			load_pdb_atoms(pdb_line,atoms, &line_atm);
		}
	}
	free(fname);
	fclose(pdbfile);
}

void load_pdb_file_without_num_atom(pdb_atom_t *atoms, pdb_seqres_t *seqres,
		const char *path_file_name){
	FILE *pdbfile=NULL;
	char pdb_line [MAX_LINE_PDB];
	char *aux = NULL;
	int line_atm = -1;
	pdbfile = open_file(path_file_name, fREAD);
	while ( fgets(pdb_line,MAX_LINE_PDB,pdbfile) != NULL){
		if (strncmp(pdb_line,"ATOM",4) == 0){
			line_atm = line_atm + 1;
			load_pdb_atoms(pdb_line,atoms, &line_atm);
		}
	}
	fclose(pdbfile);
}

int get_num_atom(const char *path_PDB_file_name){
	/*This function reads a PDB file and return the number of atoms*/
	int num_atom = 0;
	char pdb_line [MAX_LINE_PDB];
	FILE *pdbfile=NULL;
	pdbfile = open_file(path_PDB_file_name, fREAD);
	while ( fgets(pdb_line,MAX_LINE_PDB,pdbfile) != NULL){
		if (strncmp(pdb_line,"ATOM",4) == 0){
			num_atom = num_atom + 1;
		}
	}
	fclose(pdbfile);
	return num_atom;
}



