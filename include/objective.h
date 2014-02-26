#ifndef OLD_OBJECTIVE_H
#define OLD_OBJECTIVE_H

#include "enums.h"

type_fitness_energies_t str2type_objective(char *name_fitness_energies);
void type_fitness_energies2str(char *name_fitness_energies,
		const type_fitness_energies_t *type_fitness );

#endif
