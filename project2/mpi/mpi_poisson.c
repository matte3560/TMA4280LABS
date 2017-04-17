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

#include "mpi_poisson.h"

#include <stdio.h>

poisson_result_t mpi_poisson(int n)
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
	mpi_grid(grid, h, n);

	/*
	 * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
	 * defined Chapter 9. page 93 of the Lecture Notes.
	 * Note that the indexing starts from zero here, thus i+1.
	 */
	double *diag = mk_1D_array(mpi_padded_size(m), false);
	mpi_diag(diag, m, n);

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
	mpi_gen_rhs(b, grid, h, m);

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
	mpi_dst(b, m, n, false);
	mpi_transpose(bt, b, m);
	mpi_dst(bt, m, n, true);

	/*
	 * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
	 */
	mpi_solve_tu(bt, diag, m);

	/*
	 * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
	 */
	mpi_dst(bt, m, n, false);
	mpi_transpose(b, bt, m);
	mpi_dst(b, m, n, true);

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

int mpi_init(int *argc, char ***argv)
{
	int result = MPI_Init(argc, argv);
	mpi_comm = MPI_COMM_WORLD;
	MPI_Comm_rank(mpi_comm, &mpi_rank);
	MPI_Comm_size(mpi_comm, &mpi_size);

	return result;
}

int mpi_finalize()
{
	return MPI_Finalize();
}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */
void mpi_transpose(double **bt, double **b, int m)
{
	/* Indexing */
	size_t start = mpi_idx_start(m);
	size_t part = mpi_idx_part(m);
	size_t padded_size = mpi_padded_size(m);

	/* Create block datatype */
	MPI_Datatype block;
	MPI_Type_vector(part, part, padded_size, MPI_DOUBLE, &block);
	MPI_Type_create_resized(block, 0, part*sizeof(double), &block);
	MPI_Type_commit(&block);

	/* Create transposed block datatype */
	MPI_Datatype transp_block;
	MPI_Type_vector(part, 1, padded_size, MPI_DOUBLE, &transp_block);
	MPI_Type_create_hvector(part, 1, sizeof(double), transp_block, &transp_block);
	MPI_Type_create_resized(transp_block, 0, part*sizeof(double), &transp_block);
	MPI_Type_commit(&transp_block);

	/* Perform transpose */
	MPI_Alltoall(
			b[0], 1, block,
			bt[0], 1, transp_block,
			mpi_comm);

	/* Free datatypes */
	MPI_Type_free(&block);
	MPI_Type_free(&transp_block);
}

void mpi_grid(double *grid, double h, int n)
{
	for (size_t i = mpi_idx_start(n+1); i < mpi_idx_end(n+1); i++) {
		grid[i] = i * h;
	}
	MPI_Allgather(
			MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
			grid, mpi_idx_part(n+1), MPI_DOUBLE,
			mpi_comm);
}

void mpi_diag(double *diag, int m, int n)
{
	for (size_t i = mpi_idx_start(m); i < mpi_idx_end(m); i++) {
		diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
	}
	MPI_Allgather(
			MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
			diag, mpi_idx_part(m), MPI_DOUBLE,
			mpi_comm);
}

void mpi_gen_rhs(double **b, double *grid, double h, int m)
{
	/* Indexing */
	size_t start = mpi_idx_start(m);
	size_t end = mpi_idx_end(m);

	for (size_t i = start; i < end; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i-start][j] = h * h * poisson_rhs(grid[i+1], grid[j+1]);
		}
	}
}

void mpi_dst(double **b, int m, int n, bool inv)
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
		for (size_t i = 0; i < size; i++) {
			fstinv_(b[i], &n, z, &nn);
		}
	} else {
		for (size_t i = 0; i < size; i++) {
			fst_(b[i], &n, z, &nn);
		}
	}
}

void mpi_solve_tu(double **b, double *diag, int m)
{
	/* Indexing */
	size_t start = mpi_idx_start(m);
	size_t end = mpi_idx_end(m);

	for (size_t i = start; i < end; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i-start][j] = b[i-start][j] / (diag[i] + diag[j]);
		}
	}
}

double mpi_u_max(double **b, int m)
{
	/* Local bounds */
	size_t size = mpi_idx_size(m);

	/* Calculate local umax */
	double u_max = 0.0;
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < m; j++) {
			u_max = MAX(b[i][j], u_max);
		}
	}

	/* Get global umax */
	double global_u_max = 0.0;
	MPI_Allreduce(&u_max, &global_u_max, 1, MPI_DOUBLE, MPI_MAX, mpi_comm);

	return global_u_max;
}

size_t mpi_idx_start(const size_t size)
{
	/* Calculate size of partiton */
	size_t part = size / mpi_size;

	if (size % mpi_size != 0) {
		/* Size can not be evenly divided across processes */
		part++;
	}

	return part * mpi_rank;
}

size_t mpi_idx_end(const size_t size)
{
	/* Calculate size of partiton */
	size_t part = size / mpi_size;
	size_t end;

	if (size % mpi_size == 0) {
		end = part * (mpi_rank + 1);
	} else {
		/* Size can not be evenly divided across processes */
		part++;
		end = MIN(part * (mpi_rank + 1), size); // Dont go out of bounds on the highest rank process
	}

	return end;
}

size_t mpi_idx_size(const size_t size)
{
	size_t local_size;
	size_t start = mpi_idx_start(size);
	size_t end = mpi_idx_end(size);

	/* Handle edge case for small array size and large number of processes */
	if ( start > end ) {
		local_size = 0;
	} else {
		local_size = end - start;
	}

	return local_size;
}

size_t mpi_idx_part(const size_t size)
{
	size_t part = size / mpi_size;
	if ( size % mpi_size != 0 )
		part++;

	return part;
}

size_t mpi_padded_size(const size_t size)
{
	return mpi_idx_part(size) * mpi_size;
}

void mpi_allgather_mat(double** gmat, double** lmat, int m, int n)
{
	/* Indexing */
	size_t part = mpi_idx_part(m);
	size_t padded_size = mpi_padded_size(n);

	/* Distribute to all processes */
	MPI_Allgather(
			lmat[0], part * padded_size, MPI_DOUBLE,
			gmat[0], part * padded_size, MPI_DOUBLE,
			mpi_comm);
}

void mpi_allgather_vec(double* gvec, double* lvec, int size)
{
	/* Indexing */
	size_t part = mpi_idx_part(size);

	/* Distribute to all processes */
	MPI_Allgather(
			lvec, part, MPI_DOUBLE,
			gvec, part, MPI_DOUBLE,
			mpi_comm);
}

void mpi_gather_mat(double** gmat, double** lmat, int m, int n)
{
	/* Indexing */
	size_t part = mpi_idx_part(m);
	size_t padded_size = mpi_padded_size(n);

	/* Gather on root process */
	MPI_Gather(
			lmat[0], part * padded_size, MPI_DOUBLE,
			(gmat == NULL) ? NULL : gmat[0], part * padded_size, MPI_DOUBLE,
			0, mpi_comm);
}
