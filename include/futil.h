enum mode_files{fREAD,fWRITE,fAPPEND};
typedef enum mode_files mode_files_t;

#include "enums.h"

char * path_join_file(const char *path, const char *f);
/*Returns path_file which is union of path and file.
 * There is not check about path and file. Thus,
 * path must contain full path. Example:
 *            /home/myname/Execute/
 * */

FILE * open_file(__const char * filename, mode_files_t mode);
boolean_t file_is_empty(FILE *file_aux);
void set_random_file_name(char *file_name,__const char *prefix,
		__const char *sufix );
boolean_t check_exists_file(const char *path_file_name);
