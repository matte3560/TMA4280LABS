#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include "zeta.h"


int main(int argc, char** argv)
{
	int32_t n;
	const char* output = "result.txt";

	/* Check args */
	if ( argc != 2 )
	{
		printf("Invalid number of arguments\n");
		exit(1);
	}

	/* Check that n is valid */
	n = strtol( argv[1], NULL, 10 );
	if ( n < 2 ) 
	{
		printf("n must be greater than 1\n");
		exit(1);
	}

	/* Calculate pi */
	double pi = zeta_pi(n);

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
	fprintf( fp, "%.10f\n", pi );

	return 0;
}
