#include <stdio.h>
#include <stdlib.h>
#include <math.h>



int main(int argc, char *argv[]){

	float *val;
	int N_VAL; //parametro
	FILE *obj;
	int i;


	float avg,sd;


	obj = fopen( argv[1] , "rb" );
	sscanf( argv[2] , "%d" , &N_VAL );



	val = (float *) malloc( N_VAL * sizeof(float) );



	avg = 0.0;
	sd = 0.0;


	// leitura e média ao mesmo tempo
	for( i=0 ; i<N_VAL ; i++ ){
		fscanf( obj , "%f" , &val[i] );
		avg += val[i];
	}
	avg /= N_VAL;

	fclose(obj);

	//sd
	for( i=0 ; i<N_VAL ; i++ ){
		sd += ( val[i] - avg )*( val[i] - avg ) ;
	}
	sd /= N_VAL; // Calcula do desvio-padrão da POPULAÇÃO. Se os dados fossem de uma amostra, teriamos que dividir por (N_VAL-1) para calcular o desvio-padrão da amostra

	sd=sqrt(sd);


	free(val);

	printf("%f\t%f",avg,sd);


return 0;
}
