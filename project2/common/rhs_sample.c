/* Right hand side function from sample code */

double poisson_rhs(double x, double y)
{
	return 2 * (y - y*y + x - x*x);
}
