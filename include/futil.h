#ifndef OLD_FUTIL_H
#define OLD_FUTIL_H

#include <stdio.h>

#include "enums.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

enum mode_files{fREAD,fWRITE,fAPPEND};
typedef enum mode_files mode_files_t;

_2PG_CARTESIAN_EXPORT
char * path_join_file(const char *path, const char *f);
/*Returns path_file which is union of path and file.
 * There is not check about path and file. Thus,
 * path must contain full path. Example:
 *            /home/myname/Execute/
 * */

_2PG_CARTESIAN_EXPORT
FILE * open_file(__const char * filename, mode_files_t mode);
boolean_t file_is_empty(FILE *file_aux);
void set_random_file_name(char *file_name,__const char *prefix,
		__const char *sufix );
boolean_t check_exists_file(const char *path_file_name);
long get_file_size(FILE *fp);
char *get_last_line(const char *filepath);

#endif
