#include "common_poisson.h"


/*
 * This function is used for initializing the right-hand side of the equation.
 * Other functions can be defined to swtich between problem definitions.
 */

double poisson_rhs(double x, double y) {
    return 2 * (y - y*y + x - x*x);
}


/*
 * The allocation of a vectore of size n is done with just allocating an array.
 * The only thing to notice here is the use of calloc to zero the array.
 */

double *mk_1D_array(size_t n, bool zero)
{
    if (zero) {
        return (double *)calloc(n, sizeof(double));
    }
    return (double *)malloc(n * sizeof(double));
}

/*
 * The allocation of the two-dimensional array used for storing matrices is done
 * in the following way for a matrix in R^(n1*n2):
 * 1. an array of pointers is allocated, one pointer for each row,
 * 2. a 'flat' array of size n1*n2 is allocated to ensure that the memory space
 *   is contigusous,
 * 3. pointers are set for each row to the address of first element.
 */

double **mk_2D_array(size_t n1, size_t n2, bool zero)
{
    // 1
    double **ret = (double **)malloc(n1 * sizeof(double *));

    // 2
    if (zero) {
        ret[0] = (double *)calloc(n1 * n2, sizeof(double));
    }
    else {
        ret[0] = (double *)malloc(n1 * n2 * sizeof(double));
    }
    
    // 3
    for (size_t i = 1; i < n1; i++) {
        ret[i] = ret[i-1] + n2;
    }
    return ret;
}

void free_2D_array(double **array)
{
	/* Free underlying array first */
	free(array[0]);
	free(array);
}
