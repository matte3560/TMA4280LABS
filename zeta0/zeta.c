#include "zeta.h"

#include <math.h>


double zeta_term( int32_t i )
{
	return 1 / pow( i, 2 );
}

double zeta_sum( int32_t n )
{
	double sum = 0;

	for ( int32_t i = 1; i <= n; i++ )
	{
		sum += zeta_term(i);
	}

	return sum;
}

double zeta_pi( int32_t n )
{
	return sqrt( zeta_sum(n) * 6 );
}

