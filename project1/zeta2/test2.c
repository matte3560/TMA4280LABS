#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#define _USE_MATH_DEFINES
#include <math.h>

#include "zeta.h"


int main(int argc, char** argv)
{
	const char* output = "vtest.txt";

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

	for (int32_t k = 1; k <= 24; k++)
	{
		/* Calculate n = 2^k */
		int32_t n = 2 << (k-1);

		/* Calculate pi and error, output to file */
		double pi = zeta_pi(n);
		double error = fabs( M_PI - pi );
		fprintf( fp, "Error (k = %2i, n = %8i): %.10f\n", k, n, error );
	}
}

