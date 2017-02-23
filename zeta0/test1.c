#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "zeta.h"


int main(int argc, char** argv)
{
	const char* output = "utest.txt";

	/* Expected output for series with n = 3 */
	double expected_series =
		  1 / pow( 1, 2 )
		+ 1 / pow( 2, 2 )
		+ 1 / pow( 3, 2 );
	double expected_pi = sqrt( expected_series*6 );

	/* Get the computed value and calculate error */
	double computed_series = zeta_sum(3);
	double computed_pi = zeta_pi(3);
	double error_series = fabs( expected_series - computed_series );
	double error_pi = fabs( expected_pi - computed_pi );

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
	fprintf( fp, "Series (n = 3)\n");
	fprintf( fp, "Expected value: %.10f\nComputed value: %.10f\nError: %.10f\n\n",
			expected_series, computed_series, error_series );
	fprintf( fp, "Pi (n = 3)\n");
	fprintf( fp, "Expected value: %.10f\nComputed value: %.10f\nError: %.10f\n",
			expected_pi, computed_pi, error_pi );
}

