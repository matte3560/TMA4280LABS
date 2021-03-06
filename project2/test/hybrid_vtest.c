#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "hybrid_poisson.h"
#include "common_mpi_test.h"

/* Function prototypes */

int main(int argc, char** argv)
{
	mpi_init(&argc, &argv);

	MPI_RANK0( printf("Running with %i processes\n", mpi_size); );

	for (int i = 2; i <= 11; i++) {
		/* Get result */
		int n = 1 << i;
		poisson_result_t result = hybrid_poisson(n);

		if (mpi_rank == 0) {
			/* Calculate error */
			double error = vtest_max_error(result.u, result.grid, result.n);

			/* Output */
			printf("(n = %4i) Time: %7.4f\t\tMax error: %f\n", n, result.time, error);
		}

		/* Free memory */
		finalize_result(&result);
	}

	return mpi_finalize();
}
