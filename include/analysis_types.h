#ifndef OLD_ANALYSIS_TYPE_H
#define OLD_ANALYSIS_TYPE_H


/** Represents a name of file or ID with values of 
* multi-objective analysis 
*/
typedef struct sowner_file{	
	// Name of file of method or an ID that identifies a solution.
	char *file_name;
	// Assigned vobjective values which are used to apply analysis
	double *obj_values;
	// Indicates front that solution is
	int front;
	// Indicates number of solutions that are dominated by me
	int number_solutions_are_dominated;
	/* Represents global classification based on a sorting. This 
	* sorting should be front and dominated solution. Therefore, after applied 
	* this sorting is set a ranking in all solution. First solution assigned 1 and so on
	*/
	int ranking; 
}owner_file_t;




#endif