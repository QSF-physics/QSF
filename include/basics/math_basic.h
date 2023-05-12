#include "math_cheap.h"

double hermite(int order, double x)
{
	double hrm;
	if(order == 0) hrm= 1;
	else if(order == 1)
		hrm= 2 * x;
	else
		hrm= 2 * x * hermite(order - 1, x) - 2 * (order - 1) * hermite(order - 2, x);
	return hrm;
}

constexpr inline double Sign(int x) { return (double)((x > 0) - (x < 0)); }
/// @brief Rises the number base to pow power
/// @tparam T power
/// @param base input number
/// @param pow power to which to rise
/// @return
template<typename T> T constexpr Power(T base, const ind& pow)
{
	if(pow == 0) return 1;
	else if(pow == 1)
		return base;
	else if(pow == 2)
		return base * base;
	else
		return base * Power(base, pow - 1);
}

template<typename T, typename... Ts> T constexpr Max(T t, Ts... ts)
{
	if constexpr(bool(sizeof...(Ts))) return std::max(t, T(Max(ts...)));
	else
		return t;
}

template<typename T, typename... Ts> T constexpr Min(T t, Ts... ts)
{
	if constexpr(bool(sizeof...(Ts))) return std::min(t, T(Min(ts...)));
	else
		return t;
}

bool isPrime(ind num)
{
	ind p= 2;
	if(num % p == 0) return false;
	p++;
	if(num % p == 0) return false;

	ind sq_num= sqrt(num);
	while(p < sq_num)
	{
		p+= 2;
		if((num % p) == 0) return false;
	}
	return true;
}