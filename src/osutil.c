#include <stdio.h>
#include <stdlib.h>

#ifndef WIN32
#include <unistd.h>
#else
#include <process.h>
#include <direct.h>
#endif

#include <string.h>

#include "defines.h"
#include "messages.h"
#include "string_owner.h"
#include "osutil.h"
#include "futil.h"


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
	char *path_file_name;
	path_file_name = path_join_file(path, filename);
	remove(path_file_name);
	free(path_file_name);
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

