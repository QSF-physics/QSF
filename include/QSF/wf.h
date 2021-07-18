
struct KineticEnergy {
	static constexpr REP rep = REP::P;
};
struct PotentialEnergy {
	static constexpr REP rep = REP::X;
};
struct Identity {
	static constexpr REP rep = REP::NONE;
};
struct TotalEnergy
{
	static constexpr REP rep = REP::NONE;
};
struct EnergyDifference
{
	static constexpr REP rep = REP::NONE;
};
struct Symmetrize
{
	static constexpr REP rep = REP::NONE;
};
struct AntiSymmetrize
{
	static constexpr REP rep = REP::NONE;
};
struct Orthogonalize
{
	static constexpr REP rep = REP::NONE;
};
struct Normalize
{
	static constexpr REP rep = REP::NONE;
};

// template <typename ... Args>
// inline double gaussian(double at, double delta, Args...args)
// {
// 	double delta_inv = 1 / delta;
// 	double norm = 1.0 * sqrt(delta_inv * sqrt(1.0 / pi));
// 	norm = (((Args)norm) *...);
// 	return exp(-0.5 * (((args - at) * (args - at) * delta_inv * delta_inv) + ...)) * norm;
// }

template <class Hamiltonian, class GridBase, size_t Components>
struct WF : Grid <GridBase, Components>
{
	using H = Hamiltonian; //Derived type!
	using MPIStrategy = typename GridBase::MPIStrategy;
	using MPIGrids = typename GridBase::MPIGrids;

	using InducedGrid = Grid <GridBase, Components>;
	using InducedGrid::post_step;
	using InducedGrid::pos;
	using InducedGrid::local_n;
	using InducedGrid::DIM;
	using InducedGrid::psi;
	using InducedGrid::dx;
	WF(Section& settings) :InducedGrid(settings) {
		logInfo("WF init");
	}

	WF(GridBase g) :InducedGrid(g) {
		logInfo("WF init");
	}


	//TODO: if no match here pass to derived class
	template <REP R, class BO, class COMP, size_t...Is>
	inline void compute(BO& bo, COMP&& c, seq<Is...>&& s)
	{
		((bo.template storeInBuffer < Is, COMP>(1.0)), ...);
	}

	template <REP R, class BO, AXIS AX, class... Op, size_t...Is>
	inline void compute(BO& bo, AVG<AX, Op...>&&, seq<Is...>&&)
	{
		using T = AVG<AX, Op...>;

		double res = 0;
		([&]
		 {
			 res = 0;
			 for (ind i = 0; i < local_n; i++)
			 {
				 res += std::norm(psi[i]) *
					 static_cast<Hamiltonian*>(this)->
					 template call < Op >(InducedGrid::template pos<R>(i));
			 }
			 bo.template storeInBuffer < Is, T>(res);
		 }(), ...);
	}

	void addToInitialStateFromEigenstate(size_t index, size_t state, double weight);

	template <ind SRC>
	void addToInitialStateFromSource(size_t state, double weight);

	template <ind SRC, ind ...I>
	void loadFromRoutine(superpos<I... > sp);

	void initPsiByGuess(ind index);
	template <REP R >

	void addWavepacket(double at = 0.0, double width = 1.0, double px = 0.0, double py = 0.0, double pz = 0.0);

	void addSine(double px, double py = 0.0, double pz = 0.0);
	void addConstantBackground(double value = 1.0);
	struct Eigen;

};


// template <class Hamiltonian, class GridBase, size_t Components>
// struct WF < Hamiltonian, GridBase, Components, MPI::Multi> :
// 	WF < Hamiltonian, GridBase, Components, MPI::Single>
// {
// 	using base = WF < Hamiltonian, GridBase, Components, MPI::Single>;
// 	using grid = typename base::grid;
// 	using MPIGrids = typename base::MPIGrids;
// 	using MPIGrids::moreFree;
// 	using MPIGrids::lessFree;
// 	using base::absorber;
// 	// using base::absorber;
// 	using base::n; using base::nn;
// 	using base::n2; using base::m;
// 	using base::L; using base::psi;

// };


	// void updateDerivedQuantities()
	// {
	// 	nn = n * n;
	// 	n2 = n / 2;
	// 	m = pow(n, DIM);
	// 	L = (n - 1) * dx;
	// 	xmin = -0.5 * L; // xmin: lower limit -L/2 of integration domain in one direction

	// 	inv_m = 1.0 / (double)m;
	// 	sqrt_inv_m = sqrt(inv_m);
	// 	rowSize = m / n;
	// 	dVm = dV * inv_m;

	// 	//Momentum space constants
	// 	dp = 2 * pi / double(n) / dx;
	// 	pmin = -pi / dx;
	// 	pmax = abs(pmin) - dp;
	// 	// dp = 2.0 * pi * (n - 1) / double(L) / double(n);//TODO: RECHECK!
	// 	dVP = pow(dp, DIM);
	// 	sqrt_dVP = sqrt(dVP);
	// 	inv_dVP = 1.0 / dVP;
	// 	sqrt_inv_dVP = sqrt(inv_dVP);
	// 	//Same as L/(n-1)/sqrt(2*pi)
	// 	// p_scale = sqrt_inv_m * sqrt_dV * sqrt_inv_dVP;
	// 	p_scale = dV / pow(sqrt(2 * pi), DIM);
	// 	// p_scale = L/(n-1)/sqrt(2*pi);
	// 	kin_scale = POW2(dp) * 0.5;
	// 	//SCALES


	// 	//IM embedding
	// 	m_down = (n / 2 - Im::n / 2); // m_down: lower limit of ground state in the new grid
	// 	m_up = m_down + Im::n - 1;    // m_up: upper limit of ground state in the new grid
	// 	// autosetTimestep();
	// 	dt2 = 0.5 * dt;
	// 	logSETUP("Timestep dt: %g, grid spacing dx: %g, dp: %g", dt, dx, dp);
	// 	logSETUP("pmin:%g p_scale: %g", pmin, p_scale);
	// 	logSETUP("Total number of nodes m: %td, grid size L: %g", m, L);
	// 	logSETUP("n: %td, Im::n (IM): %td, Src::n: %td", n, Im::n, Src::n);
	// }