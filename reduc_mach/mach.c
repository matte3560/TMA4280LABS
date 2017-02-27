#include "mach.h"

#include <math.h>
#include <stdlib.h>
#include <mpi.h>

#include <stdio.h>

extern int mpi_size, mpi_rank;

double mach_term( double x, int32_t i )
{
	return pow( -1, i-1 ) * ( pow( x, 2*i-1 ) / ( 2*i-1 ) );
}

double mach_sum( double x, int32_t n )
{
	/* Determine partition size */
	int32_t part_size = n / mpi_size;
	if ( n % mpi_size != 0 )
		part_size++;

	/* Allocate space for vector and partition */
	double* vector = NULL;
	double* vector_part;
	if ( mpi_rank == 0 )
		vector = (double*)malloc( sizeof(double) * part_size * mpi_size );
	vector_part = (double*)malloc( sizeof(double) * part_size );

	/* Generate vector elements on rank 0 */
	if ( mpi_rank == 0 )
	{
		for ( int32_t i = 0; i < n; i++ )
		{
			vector[i] = mach_term( x, i+1 ); // Element 0 holds v_1
		}
	}

	/* Distribute vector elements */
	MPI_Scatter( vector, part_size, MPI_DOUBLE, vector_part, part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	/* Calculate partial sum */
	double partial_sum = 0;
	for ( int32_t i = 0; i < part_size; i++ )
	{
		/* Make sure last rank doesn't go out of bounds */
		if ( !( ( part_size * mpi_rank + i ) < n ) )
			break;

		partial_sum += vector_part[i];
	}

	/* Perform global reduction */
	double sum;
#ifndef RECURSIVE_DOUBLING
	/* Use MPI function */
	MPI_Allreduce( &partial_sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
	double recv_sum = 0;

	for ( int32_t d = 0; d < log2( mpi_size ); d++ )
	{
		/* Calculate the rank of the partner process for this iteration */
		int32_t mask = pow( 2, d );
		int32_t partner = mpi_rank ^ mask; // XOR

		/* Exchange sums with paired process */
		MPI_Request req;
		MPI_Isend( &partial_sum, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &req );
		MPI_Recv( &recv_sum, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		MPI_Request_free( &req );

		/* Update partial sum */
		partial_sum += recv_sum;
	}

	sum = partial_sum;
#endif

	/* Free memory */
	if ( mpi_rank == 0 )
		free( vector );
	free( vector_part );

	return sum;
}


double mach_pi( int32_t n )
{
	return 4 * ( 4*mach_sum( 1.0/5.0, n ) - mach_sum( 1.0/239.0, n ) );
}

