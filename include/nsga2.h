#ifndef OLD_NSGA2_H
#define OLD_NSGA2_H

#include "nsga2_type.h"
#include "analysis_types.h"

nsga2_t * allocate_nsga2(const int *size);
void desallocate_nsga2(nsga2_t *nsga2);
void show_nsga2(const nsga2_t *nsga2, const int *size);
void save_nsga2_v1(const nsga2_t *nsga2, const int *size);
void save_nsga2_v2(const nsga2_t *nsga2, const int *size);
void set_front(nsga2_t *nsga2, const dominancep_t *dominancep, const int *size);
void set_crowding_distance(nsga2_t *nsga2, const int *size);
void set_fitness(nsga2_t *nsga2, const int *size);
void order(nsga2_t *nsga2, const int *size);
void initialize_nsga2(nsga2_t *nsga2, const dominancep_t *dominancep, const int *size, const owner_file_t *file_names);

#endif
