/**
 * C program to solve the two-dimensional Poisson equation on
 * a unit square using one-dimensional eigenvalue decompositions
 * and fast sine transforms.
 *
 * Einar M. Rønquist
 * NTNU, October 2000
 * Revised, October 2001
 * Revised by Eivind Fonn, February 2015
 */

#include "hybrid_poisson.h"

#include <stdio.h>

poisson_result_t hybrid_poisson(int n)
{
	/* Get start time */
	MPI_Barrier(mpi_comm);
	double time_start = MPI_Wtime();

	/*
	 *  The equation is solved on a 2D structured grid and homogeneous Dirichlet
	 *  conditions are applied on the boundary:
	 *  - the number of grid points in each direction is n+1,
	 *  - the number of degrees of freedom in each direction is m = n-1,
	 *  - the mesh size is constant h = 1/n.
	 */
	int m = n - 1;
	double h = 1.0 / n;

	/*
	 * Grid points are generated with constant mesh size on both x- and y-axis.
	 */
	double *grid = mk_1D_array(mpi_padded_size(n+1), false);
	hybrid_grid(grid, h, n);

	/*
	 * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
	 * defined Chapter 9. page 93 of the Lecture Notes.
	 * Note that the indexing starts from zero here, thus i+1.
	 */
	double *diag = mk_1D_array(mpi_padded_size(m), false);
	hybrid_diag(diag, m, n);

	/*
	 * Allocate the matrices b and bt which will be used for storing value of
	 * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
	 */
	double **b = mk_2D_array(mpi_idx_part(m), mpi_padded_size(m), false);
	double **bt = mk_2D_array(mpi_idx_part(m), mpi_padded_size(m), false);

	/*
	 * Initialize the right hand side data for a given rhs function.
	 * Note that the right hand-side is set at nodes corresponding to degrees
	 * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
	 */
	hybrid_gen_rhs(b, grid, h, m);

	/*
	 * Compute \tilde G^T = S^-1 * (S * G)^T (Chapter 9. page 101 step 1)
	 * Instead of using two matrix-matrix products the Discrete Sine Transform
	 * (DST) is used.
	 * The DST code is implemented in FORTRAN in fsf.f and can be called from C.
	 * The array zz is used as storage for DST coefficients and internally for 
	 * FFT coefficients in fst_ and fstinv_.
	 * In functions fst_ and fst_inv_ coefficients are written back to the input 
	 * array (first argument) so that the initial values are overwritten.
	 */
	hybrid_dst(b, m, n, false);
	mpi_transpose(bt, b, m);
	hybrid_dst(bt, m, n, true);

	/*
	 * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
	 */
	hybrid_solve_tu(bt, diag, m);

	/*
	 * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
	 */
	hybrid_dst(bt, m, n, false);
	mpi_transpose(b, bt, m);
	hybrid_dst(b, m, n, true);

	/* Gather final result on rank 0 */
	double** u = NULL;
	if (mpi_rank == 0) {
		u = mk_2D_array(mpi_padded_size(m), mpi_padded_size(m), false);
	}
	mpi_gather_mat(u, b, m, m);

	/* Free memory */
	free(diag);
	free_2D_array(b);
	free_2D_array(bt);

	/* Calculate finish time */
	MPI_Barrier(mpi_comm);
	double time_finish = MPI_Wtime();

	/* Put result in struct */
	poisson_result_t result = {
		.time = time_finish - time_start,
		.n = n,
		.grid = grid,
		.u = u
	};

	return result;
}

void hybrid_grid(double *grid, double h, int n)
{
	/* Indexing */
	size_t start = mpi_idx_start(n+1);
	size_t end = mpi_idx_end(n+1);

#pragma omp parallel for
	for (size_t i = start; i < end; i++) {
		grid[i] = i * h;
	}
	MPI_Allgather(
			MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
			grid, mpi_idx_part(n+1), MPI_DOUBLE,
			mpi_comm);
}

void hybrid_diag(double *diag, int m, int n)
{
	/* Indexing */
	size_t start = mpi_idx_start(m);
	size_t end = mpi_idx_end(m);

#pragma omp parallel for
	for (size_t i = start; i < end; i++) {
		diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
	}
	MPI_Allgather(
			MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
			diag, mpi_idx_part(m), MPI_DOUBLE,
			mpi_comm);
}

void hybrid_gen_rhs(double **b, double *grid, double h, int m)
{
	/* Indexing */
	size_t start = mpi_idx_start(m);
	size_t end = mpi_idx_end(m);

#pragma omp parallel for
	for (size_t i = start; i < end; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i-start][j] = h * h * poisson_rhs(grid[i+1], grid[j+1]);
		}
	}
}

void hybrid_dst(double **b, int m, int n, bool inv)
{
	/*
	 * This vector will holds coefficients of the Discrete Sine Transform (DST)
	 * but also of the Fast Fourier Transform used in the FORTRAN code.
	 * The storage size is set to nn = 4 * n, look at Chapter 9. pages 98-100:
	 * - Fourier coefficients are complex so storage is used for the real part
	 *   and the imaginary part.
	 * - Fourier coefficients are defined for j = [[ - (n-1), + (n-1) ]] while 
	 *   DST coefficients are defined for j [[ 0, n-1 ]].
	 * As explained in the Lecture notes coefficients for positive j are stored
	 * first.
	 * The array is allocated once and passed as arguments to avoid doings 
	 * reallocations at each function call.
	 */
	int nn = 4*n;
	double z[nn];

	/* Local bounds */
	size_t size = mpi_idx_size(m);

	if (inv) {
#pragma omp parallel for private(nn, z)
		for (size_t i = 0; i < size; i++) {
			fstinv_(b[i], &n, z, &nn);
		}
	} else {
#pragma omp parallel for private(nn, z)
		for (size_t i = 0; i < size; i++) {
			fst_(b[i], &n, z, &nn);
		}
	}
}

void hybrid_solve_tu(double **b, double *diag, int m)
{
	/* Indexing */
	size_t start = mpi_idx_start(m);
	size_t end = mpi_idx_end(m);

#pragma omp parallel for
	for (size_t i = start; i < end; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i-start][j] = b[i-start][j] / (diag[i] + diag[j]);
		}
	}
}
