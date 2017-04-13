#include "common_test.h"

#include <stdio.h>


void test_print_mat(double** mat, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
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
