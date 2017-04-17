#include <stdio.h>
#include <stdlib.h>

#include "mpi_poisson.h"

int main(int argc, char **argv)
{
	/* Initialize MPI */
	mpi_init(&argc, &argv);

    if (argc < 2) {
		if (mpi_rank == 0) {
			printf("Usage:\n");
			printf("  %s n [export_prefix]\n\n", argv[0]);
			printf("Arguments:\n");
			printf("              n: the problem size (must be a power of 2)\n");
			printf("  export_prefix: filename prefix for exported results\n");
			printf("                 (optional, if not provided no data is exported)\n");
		}
		return mpi_finalize(); // Terminate
    }

	/* Get grid size */
    int n = atoi(argv[1]);

	poisson_result_t result = mpi_poisson(n);
	if (mpi_rank == 0) {
		printf("(P = %2i, n = %5i) Time: = %f\n", mpi_size, n, result.time);
		if (argc >= 3) {
			export_result(&result, argv[2]);
		}
	}
	finalize_result(&result);

	return mpi_finalize();
}
