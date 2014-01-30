#include <string.h>
#include <stdio.h>

#include "objective.h"
#include "messages.h"

type_fitness_energies_t str2type_fitness_energies(char *name_fitness_energies){
	if (strcmp(name_fitness_energies,"Potential") == 0){
		return fit_ener_potential;
	}else if (strcmp(name_fitness_energies,"Van_der_Waals") == 0){
		return fit_ener_edw;
	}else if (strcmp(name_fitness_energies,"Electrostatic") == 0){
		return fit_ener_ele;
	}else if (strcmp(name_fitness_energies,"Hydrophobic") == 0){
		return fit_hydrophobic;
	}else if (strcmp(name_fitness_energies,"Hydrophilic") == 0){
		return fit_hydrophilic;
	}else if (strcmp(name_fitness_energies,"Total_Area") == 0){
		return fit_total_area;
	}else if (strcmp(name_fitness_energies,"Gyrate") == 0){
		return fit_gyrate;
	}else if (strcmp(name_fitness_energies,"H_Bond") == 0){
		return fit_hbond;
	}else if (strcmp(name_fitness_energies,"H_Bond_Main") == 0){
		return fit_hbond_main;
	}else if (strcmp(name_fitness_energies,"GBSA_Solvatation") == 0){
		return fit_GBSA_Solvatation;
	}else if (strcmp(name_fitness_energies,"Stride_total") == 0){
		return fit_stride_total;
	}else if (strcmp(name_fitness_energies,"Stride_helix") == 0){
		return fit_stride_helix;
	}else if (strcmp(name_fitness_energies,"Stride_beta") == 0){
		return fit_stride_beta;
	}else{
		char msg[1024];
		sprintf(msg,"The option %s did not find at type_fitness_energies_t. Please, check it!", name_fitness_energies);
		fatal_error(msg);
	}

}
