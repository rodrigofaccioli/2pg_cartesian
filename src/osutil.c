#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

#include "defines.h"
#include "messages.h"
#include "string_owner.h"
#include "osutil.h"


int get_pid(){
	/*Returnd the process id*/
	return getpid();
}

void set_current_working_directory(char *path){
	/* Gets the current working directory
	 */
	if ( getcwd(path, MAX_PATH) == NULL){
		fatal_error("Error when try to obtain the current working directory \n");
	}
	check_path(path);
}

void delete_file(const char *path, const char *filename){
	/*Deletes file*/
	char command[MAX_PATH_FILE_NAME+30];
	strcpy(command,"rm ");
	strcat(command,path);
	strcat(command,filename);
	strcat(command, " > /dev/null 2> /dev/null");
	system(command);
}

static void check_path(char *path){
	/*This function adds end of path termination.
	 * In this moment is considered runs in POSIX plataform, only. So, the
	 * terminator is /
	 */
	int len;
	len = length_char(path);
	strcat(path,"/");
	strcat(path,"\0");
}

