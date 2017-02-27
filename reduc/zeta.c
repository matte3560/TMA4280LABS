#include "zeta.h"

#include <math.h>
#include <stdlib.h>
#include <mpi.h>

extern int mpi_size, mpi_rank;

double zeta_term( int32_t i )
{
	return 1 / pow( i, 2 );
}

result_t zeta_sum( int32_t n )
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
			vector[i] = zeta_term(i+1); // Element 0 holds v_1
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
	result_t result; // Used for storing reduction timing and global sum
	MPI_Barrier( MPI_COMM_WORLD ); // Synchronize all processes before timing start
	double time_start = MPI_Wtime();
#ifndef RECURSIVE_DOUBLING
	/* Use MPI function */
	MPI_Allreduce( &partial_sum, &result.value, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
	/* Perform recursive doubling sum */
	for ( int32_t d = 0; d < log2( mpi_size ); d++ )
	{
		/* Calculate the rank of the partner process for this iteration */
		int32_t mask = pow( 2, d );
		int32_t partner = mpi_rank ^ mask; // XOR

		/* Exchange sums with paired process */
		double recv_sum;
		MPI_Request req;
		MPI_Isend( &partial_sum, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, &req );
		MPI_Recv( &recv_sum, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
		MPI_Wait( &req, MPI_STATUS_IGNORE );

		/* Update partial sum */
		partial_sum += recv_sum;
	}

	result.value = partial_sum;
#endif
	double time_finish = MPI_Wtime();
	result.reduc_time = time_finish - time_start;

	/* Free memory */
	if ( mpi_rank == 0 )
		free( vector );
	free( vector_part );

	return result;
}

result_t zeta_pi( int32_t n )
{
	result_t result = zeta_sum(n);
	result.value = sqrt( result.value * 6 );
	return result;
}

