
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

template <typename ... Args>
inline double gaussian(double at, double delta, Args...args)
{
	double delta_inv = 1.0 / delta;
	double norm = Power(sqrt(delta_inv * sqrt(inv_pi)), sizeof...(Args));
	return exp(-0.5 * (((args - at) * (args - at) * delta_inv * delta_inv) + ...)) * norm;
}

template <class Hamiltonian, class GridBase, size_t Components>
struct WF : Grid <GridBase, Components>
{
	using H = Hamiltonian; //Derived type!
	using MPIStrategy = typename GridBase::MPIStrategy;
	using MPIGrids = typename GridBase::MPIGrids;

	using InducedGrid = Grid <GridBase, Components>;
	using InducedGrid::post_step;
	using InducedGrid::pos;
	using InducedGrid::local_m;
	using InducedGrid::DIM;
	using InducedGrid::DIMC;
	using InducedGrid::psi;
	using InducedGrid::shape;
	using InducedGrid::reverse_shape;
	using InducedGrid::dx;
	using InducedGrid::local_start;

	WF(Section& settings) :InducedGrid(settings) {
		logInfo("WF init");
	}

	WF(GridBase g) :InducedGrid(g) {
		logInfo("WF init");
	}

#pragma region Computations
	//TODO: if no match here pass to derived class
	template <REP R, class BO, class COMP>
	inline void compute(BO& bo, COMP&& c)
	{
		bo.template store <COMP>(1.0);
	}

	template <MODE M, REP R, uind ... Is>
	void evolve_internal(double delta, seq<Is...>)
	{
		ind idxs[DIMC + 1] = { 0 };
		do {
			psi[idxs[DIMC]] *=
				static_cast<Hamiltonian*>(this)->
				template expOp < M, R>(
					delta * static_cast<Hamiltonian*>(this)->template call<std::conditional_t<R == REP::P, KineticEnergy, PotentialEnergy>>(InducedGrid::template pos<R>(idxs[Is] + (Is == DIMC - 1 ? local_start : 0))...));

		} while (!((
			(idxs[Is]++, idxs[Is] < reverse_shape[Is])
			? (idxs[DIMC]++, false) : (idxs[Is] = 0, true)) && ...));
		// for (ind i = 0; i < local_n; i++)
		// {
		// 	if constexpr (DIM == 1)
		// 		psi[i] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (i + local_start));
		// 	else
		// 	{
		// 		ind readInd1 = i * InducedGrid::n;
		// 		//Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to switch x<->y for DIM>1
		// 		for (ind j = 0; j < InducedGrid::n; j++)
		// 		{
		// 			ind readInd2 = readInd1 + j;
		// 			if constexpr (DIM == 2)
		// 			{
		// 				if (R == REP::X || MPI::region)
		// 					psi[readInd2] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (i + local_start, j));
		// 				else
		// 					psi[readInd2] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (j, i + local_start));
		// 			}
		// 			else
		// 			{
		// 				ind readInd2 = readInd2 * InducedGrid::n;
		// 				for (ind k = 0; k < InducedGrid::n; k++)
		// 				{
		// 					ind readInd3 = readInd2 + k;
		// 					if (R == REP::X || MPI::region)
		// 						psi[readInd3] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (i + local_start, j, k));
		// 					else
		// 						psi[readInd3] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (j, i + local_start, k));
		// 				}
		// 			}
		// 		}
		// 	}
		// }
	}
	template <MODE M, REP R>
	void evolve(double delta)
	{
		evolve_internal<M, R>(delta, fftw_aware_indices<R>());
	}

	template <REP R, class Op, uind ... Is>
	double average(seq<Is...>)
	{
		// logInfo("%td %td %td %td", shape[0], shape[1], shape[2], sizeof...(Is));
		ind idxs[DIMC + 1]{ 0 };
		double res = 0.0;
		do {
			// logInfo("This will get added %s", typeid(Op).name());
			// if (R == REP::X)
				// logInfo("%td/%td", idxs[DIMC], local_m);
				// logInfo("[%2td %2td %2td] [%2td %2td %2td] %td %td", idxs[Is]..., reverse_shape[0], reverse_shape[1], reverse_shape[2], sizeof...(Is), idxs[DIMC]);
			// if (R == REP::X)
				// logInfo("%td %g", idxs[DIMC], std::norm(psi[idxs[DIMC]]));

			if constexpr (std::is_same_v<Op, Identity>)
				res += std::norm(psi[idxs[DIMC]]);
			else
				res += static_cast<Hamiltonian*>(this)->template call < Op >(InducedGrid::template pos<R>(idxs[Is] + (Is == DIMC - 1 ? local_start : 0))...) * std::norm(psi[idxs[DIMC]]);

		} while (!((
			(idxs[Is]++, idxs[Is] < reverse_shape[Is])
			? (idxs[DIMC]++, false) : (idxs[Is] = 0, true)) && ...));
		return res * InducedGrid::template vol<R>();
	}
	/* Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to  switch
	x<->y for DIM>1 */
	template <REP R> //FIXME: should respect propagator firstREP!
	constexpr auto fftw_aware_indices()
	{	//FIXME: | !MPI::region
		if constexpr (R == REP::P && DIM > 1)
			return switch_seq<n_seq_t<DIMC>>{};
		else return n_seq<DIMC>;
	}

	template <REP R, class BO, class... Op>
	inline void compute(BO& bo, AVG<Op...>&&)
	{
		using T = AVG<Op...>;
		bo.template store<T>(average<R, Op>(fftw_aware_indices<R>())...);
	}
#pragma endregion Computations

#pragma region Operations
	template <REP R>
	inline auto immediate(Normalize)
	{
		auto res = average<R, Identity>(fftw_aware_indices<R>());
		MPI::reduceImmediataly(&res);
		InducedGrid::multiplyArray(psi, sqrt(1.0 / res));
		// return res;
	}

	template <REP R>
	inline auto immediate(Symmetrize)
	{

	}

	template <REP R>
	inline auto immediate(AntiSymmetrize)
	{

	}
	template <REP R>
	inline auto immediate(Orthogonalize)
	{
		// if (state > 0)
		// {
		// 	int lower = 0;
		// 	([&] {
		// 		if (lower < state)
		// 			amplits[lower] = _CO_PROJ<IdentityOperator<REP::X | REP::P>>::template calc<R, opt, integral_constant<size_t, Args>>();
		// 		lower++;
		// 	 }(), ...);

		// 	MPI::reduceImmediataly(amplits, amplits_size);
		// 	for (lower = 0; lower < state; lower++)
		// 	{
		// 		cplxd* ei = wf->states[lower];
		// 		for (i = 0; i < local_m; i++)
		// 			wf[i] = wf[i] - amplits[lower] * ei[i];
		// 	}
		// }
		// auto res = average<R, Identity>(fftw_aware_indices<R>());
		// MPI::reduceImmediataly(&res);
		// InducedGrid::multiplyArray(psi, sqrt(1.0 / res));
	}

	template <REP R, class BO, class... Op>
	inline void compute(BO& bo, EARLY_OPERATION<Op...>&&)
	{
		using T = EARLY_OPERATION<Op...>;

		(immediate<R>(Op{}), ...);
		// bo.template store<T>(average<R, Op>(fftw_aware_indices<R>())...);
	}
#pragma endregion

#pragma region InitialConditions
	void setConstantValue(cxd val)
	{
		for (ind i = 0;i < local_m; i++) psi[i] = val;
	}

	template <class F, REP R, bool coords, uind ... Is>
	void add(F&& f, seq<Is...>)
	{
		// _logMPI("my local_start %td", local_start);
		ind idxs[DIMC + 1] = { 0 };
		do {
			if constexpr (coords)
				psi[idxs[DIMC]] += f(InducedGrid::template pos<R>(idxs[Is] + (Is == DIMC - 1 ? local_start : 0))...);
			else psi[idxs[DIMC]] += f((idxs[Is] + (Is == DIMC - 1 ? local_start : 0))...);

		} while (!((
			(idxs[Is]++, idxs[Is] < reverse_shape[Is])
			? (idxs[DIMC]++, false) : (idxs[Is] = 0, true)) && ...));
	}

	template <class F, REP R = REP::X>
	void addUsingCoordinateFunction(F&& f)
	{
		add<F, R, true>(std::forward<F>(f), fftw_aware_indices<R>());
	}
	template <class F, REP R = REP::X>
	void addUsingNodeFunction(F&& f)
	{
		add<F, R, false>(std::forward<F>(f), fftw_aware_indices<R>());
	}

	void addFromFile()
	{

	}

#pragma endregion InitialConditions

	void addToInitialStateFromEigenstate(size_t index, size_t state, double weight);

	template <ind SRC>
	void addToInitialStateFromSource(size_t state, double weight);

	template <ind SRC, ind ...I>
	void loadFromRoutine(superpos<I... > sp);

	// template <class F>
	// void add(F f)
	// {
	// 	operator()()
	// }

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