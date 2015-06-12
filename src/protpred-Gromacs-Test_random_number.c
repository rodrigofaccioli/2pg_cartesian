#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#include "defines.h"
#include "randomlib.h"
#include "osutil.h"

#define MANY_TIMES 10

int main(){

	int t;
	int max_int;
	max_int = 100;

	printf("float random\n");
	t = 1;
	while (t <= MANY_TIMES){
		printf("%f \n", _get_float());
		t = t + 1;
	}

	printf("int random\n");
	t = 1;
	while (t <= MANY_TIMES){
		printf("%d \n", _get_int_random_number(&max_int));		
		t = t + 1;
	}

	return 0;
}