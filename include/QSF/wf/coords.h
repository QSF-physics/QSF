struct CoordinateSystem {
	ind n;
	double dx;
};

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
	ind n2;
	ind nn;
	ind m;
	double inv_m;
	double inv_nn;
	double inv_n;

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
	double vol() {
		if constexpr (R == REP::X) return dV;
		else return dVm;
	}

	static constexpr int DIM = intDIMS(Dflag);
	static constexpr DIMS D = Dflag;
	void init()
	{
		n2 = n / 2;
		nn = n * n;
		m = Power(n, DIM);
		inv_m = 1.0 / m;
		inv_nn = 1.0 / nn;
		inv_n = 1.0 / n;
		dV = Power(dx, DIM);
		dVm = dV * inv_m;
		L = (dx * (n - 1));
		xmin = -0.5 * L;
		inv_dx = 1.0 / dx;
		inv_2dx = inv_dx / 2.0;

		dp = 2.0 * pi / double(n) / dx;
		pmin = -pi / dx;
		pmax = -pmin - dp;
		dVP = Power(dp, DIM);
		kin_scale = Power(dp, 2) * 0.5;


		logSETUP("CartesianGrid init");
		logSETUP("Total number of nodes m: %td, n: %td, grid size L: %g", m, m, L);
		logSETUP("Grid spacing dx: %g, dp: %g", dx, dp);
		logSETUP("pmin:%g kin_scale: %g", pmin, kin_scale);
		logSETUP("inv_m:%g inv_nn: %g", inv_m, inv_nn);
	}

	CartesianGrid(Section& settings)
	{

		inipp::get_value(settings, "n", n);
		inipp::get_value(settings, "dx", dx);
		init();
	}

	CartesianGrid(CoordinateSystem cs) : CoordinateSystem(cs)
	{
		init();
	}


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
		else return dp * (index >= n2 ? index - n : index); //after FFTW 0 freq is at 0 index
	}
};