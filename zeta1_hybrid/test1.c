#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <mpi.h>
#include <assert.h>

#include "zeta.h"

/* Global vars for convenience */
int mpi_size, mpi_rank;


int main(int argc, char** argv)
{
	const char* output = "vtest.txt";

	/* Initialize MPI */
	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
	MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
	assert( mpi_size % 2 == 0 ); // Make sure we have a valid number of processes

	/* Output result to textfile */
	FILE* fp = NULL;
	if ( mpi_rank == 0 )
	{
		fp = fopen( output, "w" );
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
	}

	for (int32_t k = 1; k <= 24; k++)
	{
		/* Calculate n = 2^k */
		int32_t n = 2 << (k-1);

		/* Calculate runtime and error, output to file */
		MPI_Barrier( MPI_COMM_WORLD ); // Make sure all processes are synchronized before timing starts
		double start = MPI_Wtime();
		double pi = zeta_pi(n);
		double finish = MPI_Wtime();
		double error = fabs( M_PI - pi );

		if ( mpi_rank == 0 )
			fprintf( fp, "(P = %i, n = %8i) Error: %.10f\tTime: %f\n", mpi_size, n, error, (finish - start) );
	}

	/* Terminate MPI env */
	MPI_Finalize();

	return 0;
}

