template <DIMS DIM>
struct borVec
{
	double dat[DIM];

	template <typename... Args>
	borVec(Args...args) : dat{ args... }
	{
		static_assert(DIM == sizeof...(Args), "borVec expects #DIM arguments");
		double norm = sqrt(((args * args) + ... + 0.0));
		if (norm > 0.001) for (DIMS i = 0; i < DIM;i++)
			dat[i] /= norm;
	}

	borVec() = default;

	double& operator[](int index)
	{
		return dat[index];
	}
	double operator[](int index) const
	{
		return dat[index];
	}
	borVec operator+(const borVec& a)
	{
		for (int i = 0; i < DIM; i++) dat[i] += a[i];
		return *this;
	}
	//The side of multiplications doesn't matter here

	double operator*(const borVec& a) const
	{
		double res = 0.0;
		for (int i = 0; i < DIM; i++)
			res += dat[i] * a[i];
		return res;
	}

	// borVec<DIMS> operator-() const
	// {

	// 	return { -dat[0],-dat[1],-dat[2] };

	// 	return { -dat[0],-dat[1] };

	// 	return { -dat[0] };

	// }
};
