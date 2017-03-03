#include "zeta.h"

#include <math.h>
#include <stdlib.h>
#include <mpi.h>

extern int mpi_size, mpi_rank;

double zeta_term( int32_t i )
{
	return 1 / pow( i, 2 );
}

double zeta_sum( int32_t n )
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
#pragma omp parallel for
		for ( int32_t i = 0; i < n; i++ )
		{
			vector[i] = zeta_term(i+1); // Element 0 holds v_1
		}
	}

	/* Distribute vector elements */
	MPI_Scatter( vector, part_size, MPI_DOUBLE, vector_part, part_size, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	/* Calculate partial sum */
	double partial_sum = 0;
#pragma omp parallel for reduction(+:partial_sum)
	for ( int32_t i = 0; i < part_size; i++ )
	{
		/* Make sure last rank doesn't go out of bounds */
		if ( !( ( part_size * mpi_rank + i ) < n ) )
			continue; // Cant break out of OpenMP loop

		partial_sum += vector_part[i];
	}

	/* Gather partial sums in rank 0 */
	double* vector_part_sums = NULL;
	if ( mpi_rank == 0 )
		vector_part_sums = (double*)malloc( sizeof(double) * mpi_size );
	MPI_Gather( &partial_sum, 1, MPI_DOUBLE, vector_part_sums, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD );

	/* Calculate global sum */
	double sum = 0;
	if ( mpi_rank == 0 )
	{
#pragma omp parallel for reduction(+:sum)
		for ( int32_t i = 0; i < mpi_size; i++ )
			sum += vector_part_sums[i];
	}

	/* Free memory */
	if ( mpi_rank == 0 )
	{
		free( vector );
		free( vector_part_sums );
	}
	free( vector_part );

	return sum;
}

double zeta_pi( int32_t n )
{
	return sqrt( zeta_sum(n) * 6 );
}

