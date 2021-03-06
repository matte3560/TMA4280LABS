#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "serial_poisson.h"
#include "common_test.h"

/* Function prototypes */
void test_transpose();
void test_u_max();
void test_grid();
void test_diag();
void test_solve_tu();
void test_gen_rhs();
void test_dst();
void test_dstinv();

int main(int argc, char** argv)
{
	TESTPRINT(Testing transpose function);
	test_transpose();
	TESTPASS();

	TESTPRINT(Testing u_max function);
	test_u_max();
	TESTPASS();

	TESTPRINT(Testing grid generation function);
	test_grid();
	TESTPASS();

	TESTPRINT(Testing diagonal matrix generation function);
	test_diag();
	TESTPASS();

	TESTPRINT(Test function used to calculate tilde U);
	test_solve_tu();
	TESTPASS();

	TESTPRINT(Test function used generate the right hand side);
	test_gen_rhs();
	TESTPASS();

	TESTPRINT(Test function used for discrete sine transform);
	test_dst();
	TESTPASS();

	TESTPRINT(Test function used for inverse discrete sine transform);
	test_dstinv();
	TESTPASS();

	return 0;
}


void test_transpose()
{
	/* Create simple 3x3 matrix to transpose */
	double** mat = mk_2D_array(3, 3, false);
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			mat[i][j] = 3*i + j;
		}
	}

	/* Perform transpose */
	double** mat_t = mk_2D_array(3, 3, false);
	serial_transpose(mat_t, mat, 3);

	puts("Before");
	test_print_mat(mat, 3, 3);
	puts("After");
	test_print_mat(mat_t, 3, 3);

	/* Check that the values are correct */
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			APPROX(mat_t[i][j], 3*j + i);
		}
	}

	/* Free memory */
	free_2D_array(mat);
	free_2D_array(mat_t);
}

void test_u_max()
{
	/* Create 3x3 test matrix */
	double** mat = mk_2D_array(3, 3, true);
	mat[0][0] = 1;
	mat[2][1] = 3;
	mat[2][0] = 7;
	mat[1][2] = 10;
	double u_max = serial_u_max(mat, 3);

	test_print_mat(mat, 3, 3);
	printf("u_max = %f\n", u_max);

	assert( u_max == 10 );

	/* Free memory */
	free_2D_array(mat);
}

void test_grid()
{
	int n = 4;
	double vec[n+1];
	double h = 1.0/n;
	serial_grid(vec, h, n);

	test_print_vec(vec, n+1);

	/* Check for expected values */
	APPROX(vec[0], 0.0);
	APPROX(vec[1], 1.0/4.0);
	APPROX(vec[2], 2.0/4.0);
	APPROX(vec[3], 3.0/4.0);
	APPROX(vec[4], 4.0/4.0);
}

void test_diag()
{
	int n = 4;
	int m = n-1;
	double diag[m];
	serial_diag(diag, m, n);

	test_print_vec(diag, m);

	/* Check for expected values */
	APPROX(diag[0], 2.0 * (1.0 - cos((0+1) * PI / 4)));
	APPROX(diag[1], 2.0 * (1.0 - cos((1+1) * PI / 4)));
	APPROX(diag[2], 2.0 * (1.0 - cos((2+1) * PI / 4)));
}

void test_solve_tu()
{
	/* Create matrix */
	int n = 3;
	int m = n-1;
	double** mat = mk_2D_array(m, m, false);
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			mat[i][j] = m*i + j;
		}
	}

	/* Create diagonal */
	double diag[m];
	serial_diag(diag, m, n);

	puts("Before");
	test_print_mat(mat, m, m);
	serial_solve_tu(mat, diag, m);
	puts("After");
	test_print_mat(mat, m, m);

	/* Check for correct values */
	APPROX(mat[0][0], 0);
	APPROX(mat[0][1], 1/(diag[0] + diag[1]));
	APPROX(mat[1][0], 2/(diag[1] + diag[0]));
	APPROX(mat[1][1], 3/(diag[1] + diag[1]));

	/* Free memory */
	free_2D_array(mat);
}

void test_gen_rhs()
{
	/* Create grid */
	int n = 4;
	int m = n-1;
	double vec[n+1];
	double h = 1.0/n;
	serial_grid(vec, h, n);

	/* Create matrix */
	double** mat = mk_2D_array(m, m, false);
	serial_gen_rhs(mat, vec, h, m);

	test_print_mat(mat, m, m);

	/* Check for correct values */
	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			APPROX(mat[i][j], h * h * poisson_rhs(vec[i+1], vec[j+1]));
		}
	}
}

void test_dst()
{
	/* Create 3x5 matrix to test with */
	int n = 4; // fst_ fortran function uses n-1
	int m = n-1;
	double r1[3] = { 1, 0, 1};
	double r2[3] = { 0, 2, 1};
	double r3[3] = {.5, 0, 0};
	double* mat[3] = {r1, r2, r3};
	double** mat_copy = dup_2D_array(mat, 3, 3);
	assert( mat_copy[0] != mat[0] );

	/*
	 * Perform transformation on copied matrix using code copied directly from provided sample code.
	 * This code is assumed to be correct.
	 */
	int nn = 4*n;
	double z[nn];
    for (size_t i = 0; i < m; i++) {
        fst_(mat_copy[i], &n, z, &nn);
    }

	puts("Before");
	test_print_mat(mat, m, n-1);

	serial_dst(mat, m, n, false);

	puts("After");
	test_print_mat(mat, m, n-1);

	/* Check if values are correct */
	test_cmp_mat(mat, mat_copy, m, n-1);

	/* Free memory */
	free_2D_array(mat_copy);
}

void test_dstinv()
{
	/* Create 3x5 matrix to test with */
	int n = 4; // fst_ fortran function uses n-1
	int m = n-1;
	double r1[3] = { 1, 0, 1};
	double r2[3] = { 0, 2, 1};
	double r3[3] = {.5, 0, 0};
	double* mat[3] = {r1, r2, r3};
	double** mat_copy = dup_2D_array(mat, 3, 3);
	assert( mat_copy[0] != mat[0] );

	/*
	 * Perform transformation on copied matrix using code copied directly from provided sample code.
	 * This code is assumed to be correct.
	 */
	int nn = 4*n;
	double z[nn];
    for (size_t i = 0; i < m; i++) {
        fstinv_(mat_copy[i], &n, z, &nn);
    }

	puts("Before");
	test_print_mat(mat, m, n-1);

	serial_dst(mat, m, n, true);

	puts("After");
	test_print_mat(mat, m, n-1);

	/* Check if values are correct */
	test_cmp_mat(mat, mat_copy, m, n-1);

	/* Free memory */
	free_2D_array(mat_copy);
}
