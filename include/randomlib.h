#ifndef OLD_RAMDOMLIB_H
#define OLD_RAMDOMLIB_H

long int _get_seed();
int _get_int_random_number(__const int *max_number);
double _get_double_random_number();
double _get_double_gama(const int *a, const int *b);
double _get_double_gauss(const int *g);
void _finish_random_gsl();
float _get_float();
float _get_float_max(const float *max);
float _get_float_random_interval(const float *min, const float *max);

#endif
