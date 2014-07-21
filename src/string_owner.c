#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "messages.h"
#include "defines.h"
#include "enums.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

void int2str(char *str, const int *number){
	sprintf(str,"%d",*number);
}

int str2int(__const char *str){
	return atoi(str);
}

float str2float(__const char *str){
	return atof(str);
}

double str2double(__const char *str){
	return strtod(str, NULL);
}

int length_char(__const char *str){
	return strlen(str);
}

void append_char(char *dest, __const char *source, int max_num){
	strncat(dest, source, max_num);
}

void substring(char *dest, const char __source[], int begin, int end){
/* This function copy some part from source to dest
 * The length of dest variable must be (end - begin)+1
 * This function was based on research from internet (google). However, I found
 * examples which work with printf command. Therefore, I had idea change printf
 * to sprintf command.
 */
	sprintf(dest,"%.*s", end - begin, &__source[begin]);
	strcat(dest,"\0");
}

_2PG_CARTESIAN_EXPORT
void remove_character(char *str, const char ch){
	/*This function removes ch character from str*/
	  char *tr = '\0';
	  int c;
	  int nul;

	  if (!str)
	    return;

	  tr = strdup (str);
	  c  = 0;
	  while ((tr[c] == ch))
	    c++;

	  strcpy (str,tr+c);

	  if (!str)
	    return;

	  nul = strlen(str)-1;
	  while ((nul > 0) && ((str[nul] == ch) ) ) {
	    str[nul] = '\0';
	    nul--;
	  }
	  free (tr);
}

void remove_character_enter(char *str){
	/*This function removes ch character from str*/
	  char garbage = '\n';

    char *src, *dst;
    for (src = dst = str; *src != '\0'; src++) {
        *dst = *src;
        if (*dst != garbage) dst++;
    }
    *dst = '\0';

}


void ltrim (char *str){
/*This function was based on Gromacs 4.5.3*/
  char *tr = '\0';
  int c;

  if (!str)
    return;

  tr = strdup (str);
  c  = 0;
  while ((tr[c] == ' ') || (tr[c] == '\t'))
    c++;

  strcpy (str,tr+c);
  free (tr);
}

void rtrim (char *str){
/*This function was based on Gromacs 4.5.3*/
  int nul;

  if (!str)
    return;

  nul = strlen(str)-1;
  while ((nul > 0) && ((str[nul] == ' ') || (str[nul] == '\t')) ) {
    str[nul] = '\0';
    nul--;
  }
}

_2PG_CARTESIAN_EXPORT
void trim (char *str){
/*This function was based on Gromacs 4.5.3*/
  ltrim (str);
  rtrim (str);
}

boolean_t is_equal(const char *c1, const char *c2){
	/* Compare two chars.
	 * Returns true when those chars are equal. Otherwise, returns false.
	 */
	if ( strcmp(c1,c2) ==0 ){
		return btrue;
	}else{
		return bfalse;
	}
}

/** Returns True when letter is an alphabetic character. Otherwise, returns False
*/
boolean_t is_letter(const char letter){
	if (isalpha(letter) > 0){
		return btrue;
	}
	return bfalse;
}
