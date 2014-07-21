#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
void fatal_error(const char *__restrict __message);

_2PG_CARTESIAN_EXPORT
void display_msg(const char *__restrict __message);
