#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "serial_poisson.h"
#include "mpi_poisson.h"
#include "common_mpi_test.h"

/* Function prototypes */
void test_dst(int n);

int main(int argc, char** argv)
{
	mpi_init(&argc, &argv);

	MPI_RANK0( printf("Running tests with %i processes\n", mpi_size); );

	MPI_TESTPRINT(Testing DST with small matrix);
	test_dst(8);
	MPI_TESTPASS();

	MPI_TESTPRINT(Testing DST with larger matrix);
	test_dst(512);
	MPI_TESTPASS();

	return mpi_finalize();
}


void test_dst(int n)
{
	/* Create matrix to test with */
	int m = n-1;
	double h = 1.0/n;
	double** mat = mk_2D_array(mpi_padded_size(m), mpi_padded_size(m), false);
	double grid[n+1];
	serial_grid(grid, h, n);
	serial_gen_rhs(mat, grid, h, m);

	MPI_RANK0(
		printf("n = %i, m = %i\n", n, m);
		printf("Using %ix%i matrix.\n", m, m);
		if (n<=8) {
			puts("Before");
			test_print_mat(mat, m, m);
		} else {
			puts("Matrix too big to print...");
		}
	);

	/*
	 * Create copy of matrix and apply correct serial version.
	 * (This assumes that the unit tests for the serial implementation
	 * has been successfully run first)
	 */
	double** c_mat = dup_2D_array(mat, m, m);
	serial_dst(c_mat, m, n, false);

	/* Apply parallel version and compare results */
	mpi_dst(mat, m, n, false);

	if (n<=8) {
		MPI_RANK0(
			puts("After DST");
			test_print_mat(mat, m, m);
		);
	}

	test_cmp_mat(mat, c_mat, m, m);

	/* Also test inverse DST */
	serial_dst(c_mat, m, n, true);
	mpi_dst(mat, m, n, true);

	if (n<=8) {
		MPI_RANK0(
			puts("After DST and inverse DST");
			test_print_mat(mat, m, m);
		);
	}

	test_cmp_mat(mat, c_mat, m, m);
}
