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

	double u_max = mpi_poisson(n);
	if ( mpi_rank == 0) {
		printf("u_max = %e\n", u_max);
	}

	return mpi_finalize();
}
