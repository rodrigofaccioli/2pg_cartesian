#include <string.h>
#include <stdio.h>


#include "aminoacids.h"
#include "defines.h"
#include "messages.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

static void initialize_primary_seq(primary_seq_t* seq){

	for (int i = 0; i < seq->num_res; i++){
		seq->seq_res[i].id_1 = Malloc(char,2);
		seq->seq_res[i].id_3 = Malloc(char,4);
		strcpy(seq->seq_res[i].id_1, "");
		strcpy(seq->seq_res[i].id_3, "");
		seq->seq_res[i].id = aNR;		
	}
}

primary_seq_t* allocate_primary_seq(const int *num_res){
	primary_seq_t* aux;

	aux = Malloc(primary_seq_t, 1);
	aux->num_res = *num_res;
	aux->seq_res = Malloc(amino_t, *num_res);

	initialize_primary_seq(aux);

	return aux;
}

_2PG_CARTESIAN_EXPORT
void desallocate_primary_seq(primary_seq_t* seq){
	if (seq->seq_res != NULL){
		for (int i = 0; i< seq->num_res; i++){
			free(seq->seq_res[i].id_1);
			free(seq->seq_res[i].id_3);			
		}
		free(seq->seq_res);
	}
	free(seq);
}

/** Converts string to type_terminal_charge type
* s_term_charge means the string that will be converted to type_terminal_charge
* returns type_terminal_charge
*/
type_terminal_charge_t str2terminal_charge(const char *s_term_charge){
	char *s_none, *s_ace, *s_nme;

	s_none = Malloc(char, 5); 
	s_ace  = Malloc(char, 4); 
	s_nme = Malloc(char, 4); 

	strcpy(s_none, "none");
	strcpy(s_ace, "ACE");
	strcpy(s_nme, "NME");

	type_terminal_charge_t ret;
	if ( strcmp(s_term_charge, s_none) == 0){
		ret = term_charge_none;
	}else if ( strcmp(s_term_charge,s_ace)  == 0) {
		ret = term_charge_ACE;
	}else if ( strcmp(s_term_charge, s_nme) == 0){
		ret =  term_charge_NME;
	}else{
 		char msg[300];
		snprintf(msg, sizeof(msg), "%s it is not found at  str2terminal_charge function \n", s_term_charge);
  		fatal_error(msg);		
	}

	free(s_none);
	free(s_ace);
	free(s_nme);

	return ret;
}

type_aminos_t _get_amino_id_3(char *c){
/*Receives an amino (char) and returns your id*/
		type_aminos_t amino_id;
		char msg[30];
		if ( (strcmp(c,"GLY") == 0) ){
			amino_id =  aGLY;
		}else if ( (strcmp(c,"ARG") == 0) ){
			amino_id =  aARG;
		}else if ( (strcmp(c,"ALA") == 0) ){
			amino_id =  aALA;
		}else if ( (strcmp(c,"SER") == 0) ){
			amino_id =  aSER;
		}else if ( (strcmp(c,"THR") == 0) ){
			amino_id =  aTHR;
		}else if ( (strcmp(c,"CYS") == 0) ){
			amino_id =  aCYS;
		}else if ( (strcmp(c,"VAL") == 0) ){
			amino_id =  aVAL;
		}else if ( (strcmp(c,"LEU") == 0) ){
			amino_id =  aLEU;
		}else if ( (strcmp(c,"ILE") == 0) ){
			amino_id =  aILE;
		}else if ( (strcmp(c,"MET") == 0) ){
			amino_id =  aMET;
		}else if ( (strcmp(c,"PRO") == 0) ){
			amino_id =  aPRO;
		}else if ( (strcmp(c,"PHE") == 0) ){
			amino_id =  aPHE;
		}else if ( (strcmp(c,"TYR") == 0) ){
			amino_id =  aTYR;
		}else if ( (strcmp(c,"TRP") == 0) ){
			amino_id =  aTRP;
		}else if ( (strcmp(c,"ASP") == 0) ){
			amino_id =  aASP;
		}else if ( (strcmp(c,"GLU") == 0) ){
			amino_id =  aGLU;
		}else if ( (strcmp(c,"ASN") == 0) ){
			amino_id =  aASN;
		}else if ( (strcmp(c,"GLN") == 0) ){
			amino_id =  aGLN;
		}else if ( (strcmp(c,"HIS") == 0) ){
			amino_id =  aHIS;
		}else if ( (strcmp(c,"LYS") == 0) ){
			amino_id =  aLYS;
		}else{			
			if (strcmp(c,"") == 0){
				snprintf(msg, sizeof(msg), "Amino not found, because amino variable is empty. Check it. \n");
			}else{
				snprintf(msg, sizeof(msg), "%s Amino not found. Check it. \n",c);
			}
			fatal_error(msg);
		}
		return amino_id;
}


type_aminos_t _get_amino_id_1(char c){
/*receives an amino (char) and returns its id*/
	    type_aminos_t amino_id;
		switch (c)   {
			case 'A': //alanina
				amino_id =  aALA;
				break;
			case 'V': //valina
				amino_id =  aVAL;
				break;
			case 'F': //fenilalanina
				amino_id =  aPHE;
				break;
			case 'P': //prolina
				amino_id =  aPRO;
				break;
			case 'L': //leucina
				amino_id =  aLEU;
				break;
			case 'I': //iso leucina
				amino_id =  aILE;
				break;
			case 'R': //argenina
				amino_id =  aARG;
				break;
			case 'D': //acido aspartico
				amino_id =  aASP;
				break;
			case 'E': //acido glutamico
				amino_id =  aGLU;
				break;
			case 'S': //serina
				amino_id =  aSER;
				break;
			case 'T': //treonina
				amino_id =  aTHR;
				break;
			case 'C': //cisteina
				amino_id =  aCYS;
				break;
			case 'N': //asparagina
				amino_id =  aASN;
				break;
			case 'Q': //glutanima
				amino_id =  aGLN;
				break;
			case 'H': //histidina
				amino_id =  aHIS;
				break;
			case 'K': //lisina
				amino_id =  aLYS;
				break;
			case 'Y': //tirosina
				amino_id =  aTYR;
				break;
			case 'M': //metionina
				amino_id =  aMET;
				break;
			case 'W': //triptofano
				amino_id =  aTRP;
				break;
			case 'G': //glicina
				amino_id =  aGLY;
				break;
			case 'X': //Any element
				amino_id =  aX;
				break;				
			case '0': //ACE
				amino_id =  aACE;
				break;				
			case '1': //NME
				amino_id =  aNME;
				break;				
			default:
				amino_id = aNR;
				char mens [] = "Amino value is not known. Please, check your file.";
				fatal_error(mens);
		}
		return amino_id;
}


/** set_amino_id_3 sets the name of aminoacids based on its type.
*/
void set_amino_id_3(char *amino_name, const type_aminos_t *amino_id){

		if ( *amino_id ==  aGLY ){
			strcpy(amino_name,"GLY");
		}else if ( *amino_id ==  aARG ){
			strcpy(amino_name,"ARG");
		}else if ( *amino_id ==  aALA){
			strcpy(amino_name,"ALA");
		}else if ( *amino_id ==  aSER){
			strcpy(amino_name,"SER");
		}else if ( *amino_id ==  aTHR){
			strcpy(amino_name,"THR");
		}else if ( *amino_id ==  aCYS){
			strcpy(amino_name,"CYS");
		}else if ( *amino_id ==  aVAL){
			strcpy(amino_name,"VAL");
		}else if ( *amino_id ==  aLEU){
			strcpy(amino_name,"LEU");
		}else if (*amino_id ==  aILE){
			 strcpy(amino_name,"ILE");
		}else if ( *amino_id ==  aMET ){
			strcpy(amino_name,"MET");
		}else if ( *amino_id ==  aPRO){
			strcpy(amino_name,"PRO");
		}else if ( *amino_id ==  aPHE){
			strcpy(amino_name,"PHE");
		}else if ( *amino_id ==  aTYR){
			strcpy(amino_name,"TYR");
		}else if ( *amino_id ==  aTRP){
			strcpy(amino_name,"TRP");
		}else if ( *amino_id ==  aASP){
			strcpy(amino_name,"ASP");
		}else if ( *amino_id ==  aGLU){
			strcpy(amino_name,"GLU");
		}else if ( *amino_id ==  aASN){
			strcpy(amino_name,"ASN");
		}else if ( *amino_id ==  aGLN){
			strcpy(amino_name,"GLN");
		}else if ( *amino_id ==  aHIS ){
			strcpy(amino_name,"HIS");
		}else if ( *amino_id ==  aLYS){
			strcpy(amino_name,"LYS");
		}else if ( *amino_id ==  aX){
			strcpy(amino_name,"X");
		}else if ( *amino_id ==  aACE){
			strcpy(amino_name,"ACE");
		}else if ( *amino_id ==  aNME){
			strcpy(amino_name,"NME");
		}else{
			char msg [500];
			snprintf(msg, sizeof(msg), "In set_*amino_id_3 function *amino_id not found. Check it. \n");
			fatal_error(msg);
		}
}

/** Check N-Terminal and C-Terminal of primary sequence. When there is X, it will be
* set ACE or NME.
*
* primary_seq is the primary sequence of protein
* num_res is the number of residue of protein
*/
void check_terminal_charge(char *primary_seq, const int *num_res){
	//Checking N-Terminal
	if (primary_seq[0] == 'X'){		
		primary_seq[0] = '0'; // It Means ACE
	}
	//Checking C-Terminal
	if (primary_seq[*num_res-1] == 'X'){		
		primary_seq[*num_res-1] = '1'; // It Means NME
	}
}
