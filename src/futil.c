#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include "messages.h"
#include "defines.h"
#include "futil.h"
#include "string_owner.h"
#include "randomlib.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
char * path_join_file(const char *path, const char *f){
	char *path_file;
	int len;
	len = (strlen(path)+ strlen(f)+1);
	if (len > MAX_PATH_FILE_NAME){
		fatal_error("In path_join_file function, len variable is more than MAX_PATH_FILE_NAME. Please check it. \n");
	}
	path_file = Malloc(char,len);
    strcpy(path_file,path);
    strcat(path_file,f);
    return path_file;
}

_2PG_CARTESIAN_EXPORT
FILE * open_file(__const char * filename, mode_files_t mode){
	FILE *aux_file=NULL;
	if (mode == fREAD){
		aux_file = fopen(filename, "r+");
		if (!aux_file){
			printf("ERROR: Can not open File %s\n ",filename);
			perror("Error when trying open file \n");
		}
	}else if ((mode == fWRITE)){
		aux_file = fopen(filename, "w+");
	}else if ((mode == fAPPEND)){
		aux_file = fopen(filename, "a+");
	}
	return aux_file;
}

boolean_t file_is_empty(FILE *file_aux){
	/* Checks if file is empty.
	 * Returns true when file is empty. Otherwise, false.
	 */
	//char line_aux[MAX_LINE_FILE];
	if ( feof(file_aux)){		
		return btrue;
	}else{
		return bfalse;
	}
}

void set_random_file_name(char *file_name,__const char *prefix,
		__const char *sufix ){
	/* This function builds a random name and stores it at file_name
	 * file_name must be allocated.
	 */
	int random_number, max_int;
	char aux[MAX_FILE_NAME];
	random_number = _get_int_random_number(NULL);
	int2str(aux,&random_number);
	if (prefix != NULL){
		strcpy(file_name,prefix);
		strcat(file_name,aux);
	}else{
		strcpy(file_name,aux);
	}
	if (sufix != NULL){
		strcat(file_name,sufix);
	}
}


/** Returns true if exists a file. Otherwise, returns false.
*/
boolean_t check_exists_file(const char *path_file_name){	
    struct stat buffer;
    int exist = stat(path_file_name,&buffer);
    if(exist == 0){
    	return btrue;
    }        
	return bfalse;
}

/** Get the size of FILE
*/
long get_file_size(FILE *fp){
    long fsize = 0;

    fseek(fp,0,SEEK_END);
    fsize = ftell(fp); 
    fseek(fp,0,SEEK_SET);

    return fsize;
}

/** Get the last line of file
*/
char *get_last_line(const char *filepath){
    FILE *fp;
    char buff[4096+1];
    int size,i;
    long fsize;
    if(NULL==(fp=fopen(filepath, "r"))){
        perror("file cannot open at lastline");
        return NULL;
    }
    fsize= -1L*get_file_size(fp);
    if(size=fseek(fp, MAX(fsize, -4096L), SEEK_END)){
        perror("cannot seek");
        exit(1);
    }
    size=fread(buff, sizeof(char), 4096, fp);
    fclose(fp);
    buff[size] = '\0';
    i=size-1;
    if(buff[i]=='\n'){
        buff[i] = '\0';
    }
    while(i >=0 && buff[i] != '\n')
        --i;
    ++i;
    return strdup(&buff[i]);
}
