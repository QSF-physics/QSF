template <uind size>
struct borVec
{
	double dat[size];

	template <typename... Args>
	borVec(Args...args) : dat{ args... }
	{
		static_assert(size == sizeof...(Args), "borVec expects #DIM arguments");
		double norm = sqrt(((args * args) + ... + 0.0));
		if (norm > 0.001) for (uind i = 0; i < size;i++)
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
		for (int i = 0; i < size; i++) dat[i] += a[i];
		return *this;
	}
	//The side of multiplications doesn't matter here
	double operator*(const borVec& a) const
	{
	#if DIM == 3
		return dat[0] * a[0] + dat[1] * a[1] + dat[2] * a[2];
	#elif DIM ==2
		return dat[0] * a[0] + dat[1] * a[1];
	#else 
		return dat[0] * a[0];
	#endif
	}

	borVec operator-() const
	{
	#if DIM == 3
		return { -dat[0],-dat[1],-dat[2] };
	#elif DIM == 2
		return { -dat[0],-dat[1] };
	#else
		return { -dat[0] };
	#endif
	}
};
