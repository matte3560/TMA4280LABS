#include <stdio.h>
#include <stdlib.h>

#include "serial_poisson.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage:\n");
        printf("  poisson n\n\n");
        printf("Arguments:\n");
        printf("  n: the problem size (must be a power of 2)\n");
    }

	/* Get grid size */
    int n = atoi(argv[1]);

	double u_max = serial_poisson(n);
    printf("u_max = %e\n", u_max);

	return 0;
}
