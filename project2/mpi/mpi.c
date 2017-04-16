#include <stdio.h>
#include <stdlib.h>

#include "mpi_poisson.h"

int main(int argc, char **argv)
{
	/* Initialize MPI */
	mpi_init(&argc, &argv);

    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s n\n\n", argv[0]);
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

	/* Get grid size */
    int n = atoi(argv[1]);

	poisson_result_t result = mpi_poisson(n);
	if (mpi_rank == 0) {
		printf("(P = %2i, n = %5i) Time: = %f\n", mpi_size, n, result.time);
	}
	finalize_result(&result);

	return mpi_finalize();
}
