#ifndef OLD_SPEA2_H
#define OLD_SPEA2_H

#include "spea2_type.h"
#include "analysis_types.h"

spea2_t * allocate_spea2(const int *size);
void desallocate_spea2(spea2_t *spea2);
void show_spea2(const spea2_t *spea2, const int *size);
void save_spea2_v1(const spea2_t *spea2, const int *size);
void save_spea2_v2(const spea2_t *spea2, const int *size);
void set_raw_fitness(spea2_t *spea2, const dominancep_t *dominancep, const int *size);
void set_density(spea2_t *spea2, const int *size, const int *size_external);
void set_fitness(spea2_t *spea2, const int *size);
void order(spea2_t *spea2, const int *size);
void initialize_spea2(spea2_t *spea2, const dominancep_t *dominancep, const int *size, const owner_file_t * file_names);

#endif
