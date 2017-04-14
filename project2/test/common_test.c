#include "common_test.h"

#include <stdio.h>


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
