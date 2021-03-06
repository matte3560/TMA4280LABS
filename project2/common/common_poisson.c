#include "common_poisson.h"

#include <string.h> // memcpy???
#include <stdio.h>


double *mk_1D_array(size_t n, bool zero)
{
	/*
	 * The allocation of a vectore of size n is done with just allocating an array.
	 * The only thing to notice here is the use of calloc to zero the array.
	 */
    if (zero) {
        return (double *)calloc(n, sizeof(double));
    }
    return (double *)malloc(n * sizeof(double));
}

double **mk_2D_array(size_t n1, size_t n2, bool zero)
{
	/*
	 * The allocation of the two-dimensional array used for storing matrices is done
	 * in the following way for a matrix in R^(n1*n2):
	 * 1. an array of pointers is allocated, one pointer for each row,
	 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
	 *   is contigusous,
	 * 3. pointers are set for each row to the address of first element.
	 */

    double **ret = (double **)malloc(n1 * sizeof(double *));

    if (zero) {
        ret[0] = (double *)calloc(n1 * n2, sizeof(double));
    }
    else {
        ret[0] = (double *)malloc(n1 * n2 * sizeof(double));
    }
    
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

double **dup_2D_array(double **array, size_t n1, size_t n2)
{
	/* Create new array to copy into */
    double **ret = (double **)malloc(n1 * sizeof(double *));
	ret[0] = (double *)malloc(n1 * n2 * sizeof(double));
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }

	/* Copy contents of provided array row by row */
	for (size_t i = 0; i < n1; i++)
	{
		memcpy(ret[i], array[i], n2 * sizeof(double));
	}

    return ret;
}

void free_2D_array(double **array)
{
	/* Free underlying array first */
	if (array != NULL) {
		free(array[0]);
	}
	free(array);
}

void finalize_result(poisson_result_t* result)
{
	free( result->grid );
	free_2D_array( result->u );
}

void export_result(poisson_result_t* result, const char* name)
{
	/* Filenames */
	char grid_filename[64], u_filename[64];
	sprintf(grid_filename, "%s_grid.csv", name);
	sprintf(u_filename, "%s_u.csv", name);

	printf("Exporting results to %s and %s\n", grid_filename, u_filename);

	/* Export grid data */
	FILE* grid_file = fopen(grid_filename, "w");
	fprintf(grid_file, "%.10f", result->grid[0]); // No delimiter at start of line
	for (int i = 1; i <= result->n; i++) {
		fprintf(grid_file, "%s%.10f", DELIM, result->grid[i]);
	}
	fprintf(grid_file, "\n"); // End line

	/* Export solution data (u) */
	FILE* u_file = fopen(u_filename, "w");
	for (int i = 0; i < result->n-1; i++) {
		fprintf(u_file, "%.10f", result->u[i][0]); // No delimiter at start of line
		for (int j = 1; j < result->n-1; j++) {
			fprintf(u_file, "%s%.10f", DELIM, result->u[i][j]);
		}
		fprintf(u_file, "\n");
	}
}
