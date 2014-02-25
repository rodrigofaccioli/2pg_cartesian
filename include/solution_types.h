#ifndef OLD_SOLUTION_TYPE_H
#define OLD_SOLUTION_TYPE_H

/**
 * solution_t represents a solution
 * num_obj is the number of objectives
 * obj_values is an array of objecive values
 * representation is a pointer to something struct such as protein_t
*/
typedef struct ssolution {
   int num_obj;
   double *obj_values;
   void *representation;
 }solution_t;


#endif