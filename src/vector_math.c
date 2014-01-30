/*This file contains all routines for computing
 * vectors.
 */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"vector_types.h"
 #include"defines.h"

double mod_vector(const own_vector_t *vec){
	/*Computes the module of vector*/
	return sqrt(pow(vec->x,2)+pow(vec->y,2)+pow(vec->z,2));
}

void cross_product(own_vector_t* P, const own_vector_t *va,
		const own_vector_t *vb){
	/* Computes the cross product between va and vb.
	 * These value is saved in P_res
	 * P = AxB
	*/
	P->x = (va->y*vb->z) - (vb->y*va->z);
	P->y = (va->z*vb->x) - (vb->z*va->x);
	P->z = (va->x*vb->y) - (vb->x*va->y);
}

double p_scalar(const own_vector_t *a, const own_vector_t *b, const own_vector_t *c){
	/* Computes product scalar of 3 atoms. 
	 * It is used to compute bond angle in _compute_bond_angle function at
	 * topology.c 
	*/
	own_vector_t *ab; 
	own_vector_t *bc;

	double x, y, z, p;

	ab = Malloc(own_vector_t,1);
	bc = Malloc(own_vector_t,1);

	ab->x = b->x - a->x;
	ab->y = b->y - a->y;
	ab->z = b->z - a->z;

	bc->x = c->x - b->x;
	bc->y = c->y - b->y;
	bc->z = c->z - b->z;

	x = ab->x * bc->x;
	y = ab->y * bc->y;
	z = ab->z * bc->z;

	p = x + y + z;

	free(ab);
    free(bc);

    return p;

}

void sub_vector(own_vector_t* P,const own_vector_t *va,
		const own_vector_t *vb){
	/* Computes the difference between two vectors.
	 * P = B-A
	 */
	P->x = vb->x - va->x;
	P->y = vb->y - va->y;
	P->z = vb->z - va->z;
}

long double scalar_prod(const own_vector_t* v1, const own_vector_t* v2){
	return (v1->x*v2->x) + (v1->y*v2->y) + (v1->z*v2->z);
}


long double mod_vector_long(const own_vector_long_t *vec){
	/*Computes the module of vector*/
	return sqrt(pow(vec->x,2)+pow(vec->y,2)+pow(vec->z,2));
}

void cross_product_long(own_vector_long_t* P, const own_vector_long_t *va,
		const own_vector_long_t *vb){
	/* Computes the cross product between va and vb.
	 * These value is saved in P_res
	 * P = AxB
	*/
	P->x =  (va->y*vb->z) - (vb->y*va->z);
	P->y =  (va->z*vb->x) - (vb->z*va->x);
	P->z =  (va->x*vb->y) - (vb->x*va->y);
}

void sub_vector_long(own_vector_long_t* P,const own_vector_t *va,
		const own_vector_t *vb){
	/* Computes the difference between two vectors.
	 * P = B-A
	 */
	P->x =  vb->x - va->x;
	P->y =  vb->y - va->y;
	P->z =  vb->z - va->z;
}

long double scalar_prod_long(const own_vector_long_t* v1, const own_vector_long_t* v2){
	return (v1->x*v2->x) + (v1->y*v2->y) + (v1->z*v2->z);
}

