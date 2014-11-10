#ifndef OLD_AMINOACIDS_H
#define OLD_AMINOACIDS_H

#include "aminoacids_types.h"
#include "enums.h"

#ifdef __cplusplus
extern "C"
{
#endif

primary_seq_t* allocate_primary_seq(const int *num_res);
void desallocate_primary_seq(primary_seq_t* seq);

type_aminos_t _get_amino_id_3(char *c);
type_aminos_t _get_amino_id_1(char c);
void set_amino_id_3(char *amino_name, const type_aminos_t *amino_id);
type_terminal_charge_t str2terminal_charge(const char *s_term_charge);
void check_terminal_charge(char *primary_seq, const int *num_res);

#ifdef __cplusplus
}
#endif

#endif