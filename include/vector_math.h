#ifndef OLD_VECTOR_MATH_H
#define OLD_VECTOR_MATH_H

#include"vector_types.h"

double mod_vector(const own_vector_t *vec);
double p_scalar(const own_vector_t *a, const own_vector_t *b, const own_vector_t *c);
void cross_product(own_vector_t* P, const own_vector_t *va,
		const own_vector_t *vb);
void sub_vector(own_vector_t* P,const own_vector_t *va,
		const own_vector_t *vb);
double scalar_prod(const own_vector_t* v1, const own_vector_t* v2);
long double mod_vector_long(const own_vector_long_t *vec);
void cross_product_long(own_vector_long_t* P, const own_vector_long_t *va,
		const own_vector_long_t *vb);
void sub_vector_long(own_vector_long_t* P,const own_vector_t *va,
		const own_vector_t *vb);
long double scalar_prod_long(const own_vector_long_t* v1, const own_vector_long_t* v2);
#endif
