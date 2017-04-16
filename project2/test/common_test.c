#include "common_test.h"

#include "common_poisson.h"
#include <stdio.h>

/* Function prototypes */
static double vtest_correct(double x, double y);


void test_print_mat(double** mat, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
		{
			printf("%8.4f ", mat[i][j]);
		}
		printf("\n");
	}
}

void test_print_vec(double* vec, int size)
{
	for (int i = 0; i < size-1; i++)
	{
		printf("%.4f, ", vec[i]);
	}
	printf("%.4f\n", vec[size-1]);
}

void test_cmp_vec(double* v1, double* v2, int size)
{
	for (int i = 0; i < size; i++)
	{
		APPROX(v1[i], v2[i]);
	}
}

void test_cmp_mat(double** m1, double** m2, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		test_cmp_vec(m1[i], m2[i], cols);
	}
}

double vtest_max_error(double** mat, double* grid, int n)
{
	int m = n-1;
	double error = 0.0;

	for (int i = 0; i < m; i++) {
		for (int j = 0; j < m; j++) {
			double correct = vtest_correct(grid[i+1], grid[j+1]);
			double this_error = fabs(mat[i][j] - correct);
			error = MAX(error, this_error);
		}
	}

	return error;
}

static double vtest_correct(double x, double y)
{
	return sin(PI*x)*sin(2*PI*y);
}
