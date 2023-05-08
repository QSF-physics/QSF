

inline double ftcos_32s(double x)
{
	constexpr double c1= 0.99940307;
	constexpr double c2= -0.49558072;
	constexpr double c3= 0.03679168;
	double x2;	 // The input argument squared
	x2= x * x;
	return (c1 + x2 * (c2 + c3 * x2));
}
// Cheap sine and cosine
inline double cheap_cos(double angle)
{
	// clamp to the range 0..2pi
	angle= angle - floorf(angle * inv_twopi) * twopi;
	angle= angle > 0.f ? angle : -angle;

	if(angle < pi2) return ftcos_32s(angle);
	if(angle < pi) return -ftcos_32s(pi - angle);
	if(angle < threepi2) return -ftcos_32s(angle - pi);
	return ftcos_32s(twopi - angle);
}
inline double cheap_sin(double angle) { return cheap_cos(pi2 - angle); }