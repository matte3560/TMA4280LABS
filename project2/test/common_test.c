#include "common_test.h"

#include <stdio.h>


void test_print_mat(double** mat, int size)
{
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			printf("%8.3f ", mat[i][j]);
		}
		printf("\n");
	}
}
