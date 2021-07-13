struct CoordinateSystem {};

template <DIMS Dflag, class AT, class MG = MPI::Single, class MS = MPI::Slices>
struct CartesianGrid :CoordinateSystem
{
	static_assert(std::is_base_of_v< AbsorberType, AT>,
				  "Second CartesianGrid template argument must be type derived from AbsorberType");
	static_assert(std::is_base_of_v<MPI::Strategy, MS>,
				  "Forth CartesianGrid template argument must be type derived from MPI::Strategy");
	using MPIGrids = MG;
	using MPIStrategy = MS;
	AT absorber;
	ind n;
	ind n2;
	ind nn;
	ind m;
	double inv_m;

	double dx;
	double dV;
	double dVm;
	double L;
	double xmin;
	double inv_dx;
	double inv_2dx;

	double dp;
	double pmin;
	double pmax;
	double dVP;
	double kin_scale;

	template <REP R>
	inline double vol()
	{
		if constexpr (R == REP::X) return dV;
		else return dVm;
	}
	static constexpr int DIM = intDIMS(Dflag);
	static constexpr DIMS D = Dflag;

	CartesianGrid() {}
	explicit CartesianGrid(ind n, double dx) :
		n(n), n2(n / 2), nn(n* n), m(Power(n, DIM)),
		inv_m(1.0 / m), dx(dx), dV(Power(dx, DIM)), dVm(dV* inv_m),
		L(dx* (n - 1)), xmin(-0.5 * L), inv_dx(1.0 / dx), inv_2dx(inv_dx / 2.0),
		dp(2.0 * pi / double(n) / dx), pmin(-pi / dx), pmax(-pmin - dp),
		dVP(Power(dp, DIM)), kin_scale(Power(dp, 2) * 0.5) {}


	inline cxd scalarProductAll(cxd* from, cxd* to)
	{
		cxd over = 0.;
		for (ind i = 0; i < m; i++)
			over += std::conj(to[i]) * from[i];
		return over;
	}

	template <REP R>
	double pos(ind index)
	{
		if constexpr (R == REP::X) return xmin + index * dx;
		else return index >= n2 ? index - n : index; //after FFTW 0 freq is at 0 index
	}
};