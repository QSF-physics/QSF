struct AbsorberType {};

struct NoAbsorber : AbsorberType {};

struct CAP : AbsorberType
{
	double CAPlength = 10;
	double dx;
	double eta;
	double L2_effective;
	ind nCAP;
	CAP() {}
	CAP(double CAPlength, double dx, double L) :
		CAPlength(CAPlength),
		dx(dx),
		nCAP(ceil(CAPlength / dx)),
		eta(Power(3.0 / CAPlength, 4)),
		L2_effective(L / 2 - CAPlength) {}

	void correct(ind n, double L)
	{
		logTest((n % nCAP) == 0, "n (%td) mod (nCAP (%td)) == 0", n, nCAP);
		if (n % nCAP)
		{
			ind lCAP = nCAP - 1;
			ind rCAP = nCAP + 1;
			while ((n % lCAP) && (n % rCAP)) { lCAP--; rCAP++; }
			if ((rCAP > n / 2 - 1) && (lCAP < 1))
			{
				logError("Auto enlarging and decreasing nCAP failed.");
			}
			else if (rCAP > n / 2 - 1)
			{
				logWarning("Auto enlarging nCAP failed.");
			}
			else if (lCAP < 1)
			{
				logWarning("Auto decreasing nCAP failed.");
			}
			else
			{
				if (n % rCAP == 0)
				{
					nCAP = rCAP;
					CAPlength = nCAP * dx;
					logWarning("Auto enlarging nCAP to %td, new CAPlength %g", nCAP, CAPlength);
				}
				else if (n % lCAP == 0)
				{
					nCAP = lCAP;
					CAPlength = nCAP * dx;
					logWarning("Auto decreasing nCAP to %td, new CAPlength %g", nCAP, CAPlength);
				}
			}
		}
		L2_effective = L / 2 - CAPlength;
		eta = Power(3.0 / CAPlength, 4);
		// if (isPrime(nCAP)) logWarning("nCAP is prime, and fftw hates primes");
	}
	template <typename ...Ts>
	inline double operator()(Ts ... x)
	{
		((x = (fabs(x) - L2_effective)), ...);
		if (((x <= 0) && ...)) return 1.0;
		else
		{
			((x = x > 0 ? Power(x, 4) : 0), ...);
			return exp(-(x + ...) * eta);
		}
	}
};

