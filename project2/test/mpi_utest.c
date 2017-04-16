#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "serial_poisson.h"
#include "mpi_poisson.h"
#include "common_mpi_test.h"

/* Function prototypes */
void test_dst(int n);
void test_transpose(int n);
void test_diag(int n);
void test_grid(int n);
void test_gen_rhs(int n);

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

	MPI_TESTPRINT(Testing transpose with small matrix);
	test_transpose(8);
	MPI_TESTPASS();

	MPI_TESTPRINT(Testing transpose with large matrix);
	test_transpose(512);
	MPI_TESTPASS();

	MPI_TESTPRINT(Testing diagonal matrix generation function);
	test_diag(32);
	MPI_TESTPASS();

	MPI_TESTPRINT(Testing grid generation function);
	test_grid(32);
	MPI_TESTPASS();

	MPI_TESTPRINT(Testing RHS generation function with small size);
	test_gen_rhs(8);
	MPI_TESTPASS();

	MPI_TESTPRINT(Testing RHS generation function with large size);
	test_gen_rhs(512);
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
		printf("Using %ix%i matrix\n", m, m);
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
	mpi_allgather_mat(mat, m, m);

	MPI_RANK0(
		if (n<=8) {
			puts("After DST");
			test_print_mat(mat, m, m);
		}
	);

	test_cmp_mat(mat, c_mat, m, m);

	/* Also test inverse DST */
	serial_dst(c_mat, m, n, true);
	mpi_dst(mat, m, n, true);
	mpi_allgather_mat(mat, m, m);

	MPI_RANK0(
		if (n<=8) {
			puts("After DST and inverse DST");
			test_print_mat(mat, m, m);
		}
	);

	test_cmp_mat(mat, c_mat, m, m);

	/* Free memory */
	free_2D_array(mat);
	free_2D_array(c_mat);
}

void test_transpose(int n)
{
	/* Indexing */
	size_t padded_size = mpi_padded_size(n);

	/* Create matrix to transpose */
	double** mat = mk_2D_array(padded_size, padded_size, false);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			mat[i][j] = n*i + j;
		}
	}

	MPI_RANK0(
		printf("Using %ix%i matrix\n", n, n);
		if (n<=8) {
			puts("Before");
			test_print_mat(mat, n, n);
		} else {
			puts("Matrix too big to print...");
		}
	);

	/*
	 * Transpose using serial version to generate correct result.
	 * (This assumes that the unit tests for the serial implementation
	 * has been successfully run first)
	 */
	double** c_tmat = mk_2D_array(n, n, false);
	serial_transpose(c_tmat, mat, n);

	/* Transpose with MPI implementation and compare result */
	double** tmat = mk_2D_array(padded_size, padded_size, false);
	mpi_transpose(tmat, mat, n);

	MPI_RANK0(
		if (n<=8) {
			puts("After");
			test_print_mat(tmat, n, n);
		}
	);

	test_cmp_mat(tmat, c_tmat, n, n);

	/* Free memory */
	free_2D_array(mat);
	free_2D_array(tmat);
	free_2D_array(c_tmat);
}

void test_diag(int n)
{
	int m = n-1;

	MPI_RANK0( printf("Using n = %i, m = %i\n", n, m); );

	double diag[mpi_padded_size(m)];
	double c_diag[m];

	serial_diag(c_diag, m, n);
	mpi_diag(diag, m, n);

	MPI_RANK0(
		puts("Result");
		test_print_vec(diag, m);
	);
	
	test_cmp_vec(diag, c_diag, m);
}

void test_grid(int n)
{
	double h = 1.0 / n;

	MPI_RANK0( printf("Using n = %i, h = %f\n", n, h); );

	double grid[mpi_padded_size(n+1)];
	double c_grid[n+1];

	serial_grid(c_grid, h, n);
	mpi_grid(grid, h, n);

	MPI_RANK0(
		puts("Result");
		test_print_vec(grid, n+1);
	);

	test_cmp_vec(grid, c_grid, n+1);
}

void test_gen_rhs(int n)
{
	int m = n-1;
	double h = 1.0/n;

	MPI_RANK0( printf("Using n = %i, m = %i, h = %f\n", n, m, h); );

	/* Generate grid */
	double grid[n+1];
	serial_grid(grid, h, n);

	/* Allocate arrays */
	size_t padded_size = mpi_padded_size(m);
	double** mat = mk_2D_array(padded_size, padded_size, false);
	double** c_mat = mk_2D_array(m, m, false);

	/* Generate rhs with serial and MPI func */
	serial_gen_rhs(c_mat, grid, h, m);
	mpi_gen_rhs(mat, grid, h, m);
	mpi_allgather_mat(mat, m, m);

	MPI_RANK0(
		if (n<=8) {
			puts("Result");
			test_print_mat(mat, m, m);
		} else {
			puts("Too large to print...");
		}
	);

	/* Check for correct result */
	test_cmp_mat(mat, c_mat, m, m);

	/* Free mem */
	free_2D_array(mat);
	free_2D_array(c_mat);
}
