struct AbsorberType {};

struct NoAbsorber : AbsorberType {};

template <class GridType>
struct CAP;

template <DIMS size, class MPIRegions>
struct CAP<CartesianGrid_<size, MPIRegions>> : AbsorberType, CartesianGrid_<size, MPIRegions>
{
	using Base = CartesianGrid_<size, MPIRegions>;
	using Base::n2;

	ind nCAP;
	double eta;

	CAP(Base base, ind nCAP) : Base(base), nCAP(nCAP), eta(Power(3.0 / double(nCAP), 4)) {}
	CAP(Section& settings) : Base(settings)
	{
		inipp::get_value(settings, "nCAP", nCAP);
		eta = Power(3.0 / double(nCAP), 4);
	}

	//Takes negative distance from the corresponding edge
	template <uind ... dirs, typename ...Nodes>
	inline double absorb(double delta, Nodes ... nodes)
	{
		([&] {nodes = nCAP + ind(nodes);}(), ...);
		if (((nodes < 0) && ...)) return 1.0;
		else return exp(-delta * double(((nodes > 0 ? Power(nodes, 4) : 0.0) + ...)) * eta);
	}

	template <uind ... dirs, typename ...Nodes>
	inline double mask(Nodes ... nodes)
	{
		([&] {nodes = nCAP + ind(nodes);}(), ...);
		if (((nodes < 0) && ...)) return 1.0;
		else return exp(-double(((nodes > 0 ? Power(nodes, 4) : 0.0) + ... + 0)) * eta);
	}

	template <uind ... dirs, typename ...Nodes>
	inline double inv_mask(Nodes ... nodes)
	{
		return 1.0 - mask<dirs...>(nodes...);
	}

};

// void correct(ind n, double L)
// 	{
// 		logTest((n % nCAP) == 0, "n (%td) mod (nCAP (%td)) == 0", n, nCAP);
// 		if (n % nCAP)
// 		{
// 			ind lCAP = nCAP - 1;
// 			ind rCAP = nCAP + 1;
// 			while ((n % lCAP) && (n % rCAP)) { lCAP--; rCAP++; }
// 			if ((rCAP > n / 2 - 1) && (lCAP < 1))
// 			{
// 				logError("Auto enlarging and decreasing nCAP failed.");
// 			}
// 			else if (rCAP > n / 2 - 1)
// 			{
// 				logWarning("Auto enlarging nCAP failed.");
// 			}
// 			else if (lCAP < 1)
// 			{
// 				logWarning("Auto decreasing nCAP failed.");
// 			}
// 			else
// 			{
// 				if (n % rCAP == 0)
// 				{
// 					nCAP = rCAP;
// 					CAPlength = nCAP * dx;
// 					logWarning("Auto enlarging nCAP to %td, new CAPlength %g", nCAP, CAPlength);
// 				}
// 				else if (n % lCAP == 0)
// 				{
// 					nCAP = lCAP;
// 					CAPlength = nCAP * dx;
// 					logWarning("Auto decreasing nCAP to %td, new CAPlength %g", nCAP, CAPlength);
// 				}
// 			}
// 		}
// 		L2_effective = L / 2 - CAPlength;
// 		eta = Power(3.0 / CAPlength, 4);
// 		// if (isPrime(nCAP)) logWarning("nCAP is prime, and fftw hates primes");
// 	}