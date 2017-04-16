#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "mpi_poisson.h"
#include "common_mpi_test.h"

/* Function prototypes */
void test_index();

int main(int argc, char** argv)
{
	mpi_init(&argc, &argv);

	MPI_RANK0( printf("Running tests with %i processes\n", mpi_size); );

	MPI_TESTPRINT(Testing indexing functions);
	test_index();
	MPI_TESTPASS();

	return mpi_finalize();
}


void test_index() {
	const int m = 15;
	MPI_RANK0( printf("Using desired array size m = %i\n", m); );


	/* Check partition size */
	size_t part = mpi_idx_part(m);
	MPI_RANK0( printf("Padded partition size = %zu\n", part); );

	size_t c_part = (m % mpi_size) ? (m / mpi_size + 1) : (m / mpi_size);
	assert( part == c_part );


	/* Check padded allocation size */
	size_t padded_size = mpi_padded_size(m);
	MPI_RANK0( printf("Padded array size = %zu\n", padded_size); );

	size_t c_padded_size = part * mpi_size;
	assert( padded_size == c_padded_size );


	/* Check start index for local partition */
	MPI_RANK0( puts("Checking start index for all processes, print may be scrambled"); );
	size_t start = mpi_idx_start(m);
	MPI_SYNC( printf("Rank %i start index = %zu\n", mpi_rank, start); );

	size_t c_start = part * mpi_rank;
	assert( start == c_start );


	/* Check unpadded partition end index */
	MPI_RANK0( puts("Checking unpadded end index for all processes, print may be scrambled"); );
	size_t end = mpi_idx_end(m);
	MPI_SYNC( printf("Rank %i unpadded end index = %zu\n", mpi_rank, end); );

	size_t c_end = MIN( part * (mpi_rank + 1), m );
	assert( end == c_end );


	/* Check unpadded partition size */
	MPI_RANK0( puts("Checking unpadded partition size for all processes, print may be scrambled"); );
	size_t size = mpi_idx_size(m);
	MPI_SYNC( printf("Rank %i unpadded partition size = %zu\n", mpi_rank, size); );

	size_t c_size = start < end ? end - start : 0;
	assert( size == c_size );
}
