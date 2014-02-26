#include <string.h>
#include <stdio.h>

#include "objective.h"
#include "messages.h"

type_fitness_energies_t str2type_objective(char *name_fitness_energies){
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

void type_fitness_energies2str(char *name_fitness_energies,
		const type_fitness_energies_t *type_fitness ){
	if (*type_fitness == fit_ener_potential){
		strcpy(name_fitness_energies,"Potential");
	}else if (*type_fitness == fit_ener_edw){
		strcpy(name_fitness_energies,"Van_der_Waals");
	}else if (*type_fitness == fit_ener_ele){
		strcpy(name_fitness_energies,"Electrostatic");
	}else if (*type_fitness == fit_hydrophobic){
		strcpy(name_fitness_energies,"Hydrophobic");
	}else if (*type_fitness == fit_hydrophilic){
		strcpy(name_fitness_energies,"Hydrophilic");
	}else if (*type_fitness == fit_total_area){
		strcpy(name_fitness_energies,"Total_Area");
	}else if (*type_fitness == fit_gyrate){
		strcpy(name_fitness_energies,"Gyrate");
	}else if (*type_fitness == fit_hbond){
		strcpy(name_fitness_energies,"H_Bond");
	}else if (*type_fitness == fit_hbond_main){
		strcpy(name_fitness_energies,"H_Bond_Main");
	}else if (*type_fitness == fit_GBSA_Solvatation){
		strcpy(name_fitness_energies,"GBSA_Solvatation");
	}else if (*type_fitness == fit_stride_total){
		strcpy(name_fitness_energies,"Stride_total");
	}else if (*type_fitness == fit_stride_helix){
		strcpy(name_fitness_energies,"Stride_helix");
	}else if (*type_fitness == fit_stride_beta){
		strcpy(name_fitness_energies,"Stride_beta");
	}else{
		char msg[1024];
		sprintf(msg,"The option did not find at type_fitness_energies2str function. Please, check it!");
		fatal_error(msg);
	}
}
