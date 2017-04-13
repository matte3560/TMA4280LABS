#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "serial_poisson.h"
#include "common_test.h"

/* Function prototypes */
void test_transpose();

int main(int argc, char** argv)
{
	TESTPRINT(Testing transpose function);
	test_transpose();
}


void test_transpose()
{
	/* Create simple 3x3 matrix to transpose */
	double** mat = mk_2D_array(3, 3, false);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
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
}
