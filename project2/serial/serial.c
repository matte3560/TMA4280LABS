#include <stdio.h>
#include <stdlib.h>

#include "serial_poisson.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s n\n\n", argv[0]);
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

	/* Get grid size */
    int n = atoi(argv[1]);

	poisson_result_t result = serial_poisson(n);
	printf("(n = %i) Time: %f\n", n, result.time);
	finalize_result(&result);

	return 0;
}
