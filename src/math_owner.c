#include <math.h>

#include "defines.h"
#include "messages.h"

#include "2pg_cartesian_export.h"

_2PG_CARTESIAN_EXPORT
float degree2radians(const float *degree){
	if (*degree > 360){
		fatal_error("Invalid value for degree");
	}
	return *degree*PI/180;
}

_2PG_CARTESIAN_EXPORT
float radians2degree(const float *radians){
	return *radians*180/PI;
}

_2PG_CARTESIAN_EXPORT
float radians2degree_double(const double *radians){
	return *radians*180/PI;
}
