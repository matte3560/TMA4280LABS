#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <mpi.h>
#include <assert.h>

#include "zeta.h"

/* Global vars for convenience */
int mpi_size, mpi_rank;


int main(int argc, char** argv)
{
	int32_t n;
	const char* output = "result.txt";

	/* Initialize MPI */
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
	assert( mpi_size % 2 == 0 ); // Make sure we have a valid number of processes

	/* Check args */
	if ( argc != 2 )
	{
		printf("Invalid number of arguments\n");
		MPI_Finalize();
		exit(1);
	}

	/* Check that n is valid */
	n = strtol( argv[1], NULL, 10 );
	if ( n < 2 ) 
	{
		printf("n must be greater than 1\n");
		MPI_Finalize();
		exit(1);
	}

	/* Calculate pi */
	double pi = zeta_pi(n);

	if ( mpi_rank == 0 )
	{
		/* Output result to textfile */
		FILE* fp = fopen( output, "w" );
		if ( fp == NULL )
		{
			printf("Failed to write result to textfile\n");
			MPI_Finalize();
			exit(1);
		}
		else
		{
			printf("Writing result to %s\n", output);
		}
		fprintf( fp, "(P = %i, n = %i): %.10f\n", mpi_size, n, pi );
	}

	/* Terminate MPI env */
	MPI_Finalize();

	return 0;
}
