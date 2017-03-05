#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <omp.h>
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

		/* Calculate pi, error and runtime, output to file */
		double t_start = omp_get_wtime();
		double pi = zeta_pi(n);
		double t_finish = omp_get_wtime();
		double error = fabs( M_PI - pi );
		fprintf( fp, "(k = %2i, n = %8i) Error: %.10f\tTime: %f\n", k, n, error, (t_finish - t_start) );
	}
}

