#include <stdio.h>
#include <stdlib.h>

#include "serial_poisson.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        printf("Usage:\n");
        printf("  %s n [export_prefix]\n\n", argv[0]);
        printf("Arguments:\n");
        printf("              n: the problem size (must be a power of 2)\n");
        printf("  export_prefix: filename prefix for exported results\n");
		printf("                 (optional, if not provided no data is exported)\n");
		return 1; // Terminate
    }

	/* Get grid size */
    int n = atoi(argv[1]);

	poisson_result_t result = serial_poisson(n);
	printf("(n = %i) Time: %f\n", n, result.time);
	if (argc >= 3) {
		export_result(&result, argv[2]);
	}
	finalize_result(&result);

	return 0;
}
