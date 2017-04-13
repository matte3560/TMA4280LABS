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

#include "serial_poisson.h"

double serial_poisson(int n)
{
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
	double *grid = mk_1D_array(n+1, false);
	serial_grid(grid, h, n);

	/*
	 * The diagonal of the eigenvalue matrix of T is set with the eigenvalues
	 * defined Chapter 9. page 93 of the Lecture Notes.
	 * Note that the indexing starts from zero here, thus i+1.
	 */
	double *diag = mk_1D_array(m, false);
	serial_diag(diag, n, m);

	/*
	 * Allocate the matrices b and bt which will be used for storing value of
	 * G, \tilde G^T, \tilde U^T, U as described in Chapter 9. page 101.
	 */
	double **b = mk_2D_array(m, m, false);
	double **bt = mk_2D_array(m, m, false);

	/*
	 * Initialize the right hand side data for a given rhs function.
	 * Note that the right hand-side is set at nodes corresponding to degrees
	 * of freedom, so it excludes the boundary (bug fixed by petterjf 2017).
	 */
	serial_gen_rhs(b, grid, h, m);

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
	serial_dst(b, n, m, false);
	serial_transpose(bt, b, m);
	serial_dst(bt, n, m, true);

	/*
	 * Solve Lambda * \tilde U = \tilde G (Chapter 9. page 101 step 2)
	 */
	serial_solve_tu(bt, diag, m);

	/*
	 * Compute U = S^-1 * (S * Utilde^T) (Chapter 9. page 101 step 3)
	 */
	serial_dst(bt, n, m, false);
	serial_transpose(b, bt, m);
	serial_dst(b, n, m, true);

	/*
	 * Compute maximal value of solution for convergence analysis in L_\infty
	 * norm.
	 */
	double u_max = serial_u_max(b, m);

	/* Free memory */
	free(diag);
	free(grid);
	free_2D_array(b);
	free_2D_array(bt);

	return u_max;
}

/*
 * Write the transpose of b a matrix of R^(m*m) in bt.
 * In parallel the function MPI_Alltoallv is used to map directly the entries
 * stored in the array to the block structure, using displacement arrays.
 */

void serial_transpose(double **bt, double **b, int m)
{
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++) {
			bt[i][j] = b[j][i];
		}
	}
}

void serial_grid(double *grid, double h, int n)
{
	for (size_t i = 0; i < n+1; i++) {
		grid[i] = i * h;
	}
}

void serial_diag(double *diag, int n, int m)
{
	for (size_t i = 0; i < m; i++) {
		diag[i] = 2.0 * (1.0 - cos((i+1) * PI / n));
	}
}

void serial_gen_rhs(double **b, double *grid, double h, int m)
{
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i][j] = h * h * poisson_rhs(grid[i+1], grid[j+1]);
		}
	}
}

void serial_dst(double **b, int n, int m, bool inv)
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

	if (inv) {
		for (size_t i = 0; i < m; i++) {
			fst_(b[i], &n, z, &nn);
		}
	} else {
		for (size_t i = 0; i < m; i++) {
			fstinv_(b[i], &n, z, &nn);
		}
	}
}

void serial_solve_tu(double **b, double *diag, int m)
{
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++) {
			b[i][j] = b[i][j] / (diag[i] + diag[j]);
		}
	}
}

double serial_u_max(double **b, int m)
{
	double u_max = 0.0;
	for (size_t i = 0; i < m; i++) {
		for (size_t j = 0; j < m; j++) {
			u_max = MAX(b[i][j], u_max);
		}
	}
	return u_max;
}
