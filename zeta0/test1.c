#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "zeta.h"


int main(int argc, char** argv)
{
	const char* output = "utest.txt";

	/* Correct output for series with n = 3 */
	double correct =
		  1 / pow( 1, 2 )
		+ 1 / pow( 2, 2 )
		+ 1 / pow( 3, 2 );

	/* Get the computed value and calculate error */
	double computed = zeta_sum(3);
	double error = fabs( correct - computed );

	/* Output result to textfile */
	FILE* fp = fopen( output, "w" );
	if ( fp == NULL )
	{
		printf("Failed to write result to textfile\n");
		exit(1);
	}
	else
	{
		printf("Writing result to %s\n", output);
	}
	fprintf( fp, "Correct value: %.10f\nComputed value: %.10f\nError: %.10f\n", correct, computed, error );
}

