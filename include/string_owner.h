#include "enums.h"

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

void int2str(char *str, const int *number);
int str2int(__const char *str);
float str2float(__const char *__str);
double str2double(__const char *str);
int length_char(const char *__str);
void append_char(char *dest, __const char *source, int max_num);
void substring(char *dest, const char __source[], int begin, int end);
_2PG_CARTESIAN_EXPORT
void trim (char *str);
void rtrim (char *str);
void ltrim (char *str);
_2PG_CARTESIAN_EXPORT
void remove_character(char *str, const char ch);
void remove_character_enter(char *str);
boolean_t is_equal(const char *c1, const char *c2);
boolean_t is_letter(const char letter);
