#include <math.h>

#include "defines.h"
#include "messages.h"

float degree2radians(const float *degree){
	if (*degree > 360){
		fatal_error("Invalid value for degree");
	}
	return *degree*PI/180;
}

float degree2radians_no_pointer(float degree){
	return degree2radians(&degree);
}

float radians2degree(const float *radians){
	return *radians*180/PI;
}

float radians2degree_double(const double *radians){
	return *radians*180/PI;
}
