#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "serial_poisson.h"
#include "common_test.h"

/* Function prototypes */
void test_transpose();
void test_u_max();
void test_grid();

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
	test_print_mat(mat, 3);
	puts("After");
	test_print_mat(mat_t, 3);

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

	test_print_mat(mat, 3);
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
