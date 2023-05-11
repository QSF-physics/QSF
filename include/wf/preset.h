template <typename ... Args>
inline double gaussian(double at, double delta, Args...args)
{
	double delta_inv = 1.0 / delta;
	double norm = Power(sqrt(delta_inv * sqrt(inv_pi)), sizeof...(Args));
	return exp(-0.5 * (((args - at) * (args - at) * delta_inv * delta_inv) + ...)) * norm;
}