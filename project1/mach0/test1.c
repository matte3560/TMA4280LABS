#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "mach.h"


int main(int argc, char** argv)
{
	const char* output = "utest.txt";

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

	double x;
	/* Check series with x = 1/5 and n = 3 */
	x = 1.0/5.0;
	double expected_series1 =
		   pow( x, 2 - 1 ) / (2 - 1)
		 - pow( x, 4 - 1 ) / (4 - 1)
		 + pow( x, 6 - 1 ) / (6 - 1);

	double computed_series1 = mach_sum( x, 3 );
	double error_series1 = fabs( expected_series1 - computed_series1 );

	/* Check series with x = 1/239 and n = 3 */
	x = 1.0/239.0;
	double expected_series2 =
		   pow( x, 2 - 1 ) / (2 - 1)
		 - pow( x, 4 - 1 ) / (4 - 1)
		 + pow( x, 6 - 1 ) / (6 - 1);

	double computed_series2 = mach_sum( x, 3 );
	double error_series2 = fabs( expected_series2 - computed_series2 );

	/* Check computation of pi with n = 3 */
	double expected_pi = 4 * ( 4*computed_series1 - computed_series2 );
	double computed_pi = mach_pi( 3 );
	double error_pi = fabs( expected_pi - computed_pi );

	fprintf( fp, "Series 1 (x = 1/5, n = 3)\n");
	fprintf( fp, "Expected value: %.10f\nComputed value: %.10f\nError: %.10f\n\n",
			expected_series1, computed_series1, error_series1 );

	fprintf( fp, "Series 2 (x = 1/239, n = 3)\n");
	fprintf( fp, "Expected value: %.10f\nComputed value: %.10f\nError: %.10f\n\n",
			expected_series2, computed_series2, error_series2 );

	fprintf( fp, "Pi (n = 3)\n");
	fprintf( fp, "Expected value: %.10f\nComputed value: %.10f\nError: %.10f\n",
			expected_pi, computed_pi, error_pi );
}

