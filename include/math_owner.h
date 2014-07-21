#ifndef OLD_MATH_OWNER_H
#define OLD_MATH_OWNER_H

#ifdef WIN32
#include "2pg_cartesian_export.h"
#else
#include "2pg_cartesian_export_linux.h"
#endif

_2PG_CARTESIAN_EXPORT
float degree2radians(const float *degree);
_2PG_CARTESIAN_EXPORT
float radians2degree(const float *radians);
_2PG_CARTESIAN_EXPORT
float radians2degree_double(const double *radians);

#endif
