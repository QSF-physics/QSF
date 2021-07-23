
struct KineticEnergy {
	static constexpr REP rep = REP::P;
	static constexpr bool late = false;
};
struct PotentialEnergy {
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
};

struct Identity {
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = false;
};

struct Symmetrize
{
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = false;
};
struct AntiSymmetrize
{
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = false;
};
struct Orthogonalize
{
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = false;
};
struct Normalize
{
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = false;
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
	using InducedGrid::m_l;
	using InducedGrid::DIM;
	using InducedGrid::DIMC;
	using InducedGrid::psi;
	using InducedGrid::shape_l;
	using InducedGrid::shape_t;
	using InducedGrid::dx;
	using InducedGrid::n0_o;
	using InducedGrid::n1_o;

	WF(Section& settings) :InducedGrid(settings) {
		logInfo("WF init");
	}

	WF(GridBase g) :InducedGrid(g) {
		logInfo("WF init");
	}

#pragma region Computations
	//TODO: if no match here pass to derived class
	template <MODE M, REP R, class BO, class COMP>
	inline void compute(BO& bo, COMP&& c)
	{
		bo.template store <M, COMP>(1.0);
	}



	// ind constexpr offset<R>(ind Is)
	// {
	// 	return (Is == 0 ? n0_o : 0);
	// }
	/* Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to
	switch x<->y for DIM>1 if working on opossite rep */
	// template <ind Is>
	// ind static constexpr swap01 = Is == 0 ? 1 : (Is == 1 ? 0 : Is);

	template <ind Is>
	ind static constexpr reorder = DIMC - 1 - Is;//((R == REP::X) ? Is : swap01<Is>);
	/* When using mpi the first coordinate needs to be shifted
	NOTE: This is so when passing to functions accepting x,y,z...,
	internally the condition is different Is==0 for REP::X and Is==1 for REP::P*/
	template <REP R>
	ind constexpr offset(ind Is)
	{
		if constexpr (REP::X == R) return Is == 0 ? n0_o : 0;
		else return Is == 0 ? n1_o : 0;
	};
// {
// 	if constexpr (R == REP::X) return DIMC - 1 - Is;
// 	else return DIMC - 1 - (Is == 0 ? 1 : (Is == 1 ? 0 : Is));
// }
// static constexpr auto transposed_indices = swap_seq<n_seq_t<DIMC>>{};
// //flip_seq_t<DIMC - 1, n_seq_t<DIMC>>;
	static constexpr auto indices = n_seq<DIMC>;


	template <MODE M, REP R, uind ... Is>
	inline void evolve_(double delta, seq<Is...>)
	{
		ind counters[DIMC + 1] = { 0 };
		do {
			psi[counters[DIMC]] *= static_cast<Hamiltonian*>(this)->template expOp < M, R>
				(delta * static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(counters[Is] + offset<R>(Is))...));



			// if (R == REP::P)
			// {
			// 	psi[counters[DIMC]] = { MPI::rID,
			// 	   static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(counters[Is] + offset<R>(Is))...) };
			// }

			// if (R == REP::P)
				// _logMPI("%2d %2td [%2td %2td] (%2td %2td) [%2td %2td] %10g %10g %10g", MPI::rID, MPI::pID * m_l + counters[DIMC],
						// (counters[Is] + offset<R>(Is))..., (counters[Is])...,
						// shape_l[Is]..., std::norm(psi[counters[DIMC]]),
						// static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(counters[Is] + offset<R>(Is))...), static_cast<Hamiltonian*>(this)->template expOp < M, R>
						// (delta * static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(counters[Is] + offset<R>(Is))...)));

			// static_cast<Hamiltonian*>(this)->
				// template expOp < M, R>(
					// delta * static_cast<Hamiltonian*>(this)->template call<std::conditional_t<R == REP::P, KineticEnergy, PotentialEnergy>>(InducedGrid::template pos<R>(counters[Is] + offset<R>(Is))...));
		} while (!(...&& ((counters[reorder<Is>]++, counters[reorder<Is>] < (REP::X == R ? shape_l[reorder<Is>] : shape_t[reorder<Is>]))
						  ? (counters[DIMC]++, false) : (counters[reorder<Is>] = 0, true))));
	}
	template <MODE M, REP R>
	inline void evolve(double delta)
	{
		evolve_<M, R>(delta, indices);
	}


	template <REP R, class Op, uind ... Is>
	double average_(seq<Is...>)
	{
		// logInfo("%td %td %td %td", shape_l[0], shape_l[1], shape_l[2], sizeof...(Is));
		ind counters[DIMC + 1]{ 0 };
		double res = 0.0;
		do {
			// logInfo("This will get added %s", typeid(Op).name());
			// if (R == REP::P)
				// logInfo("%2td %2td [%2td %2td] [%2td %2td] %td %td", Is..., counters[Is]..., shape_l[Is]..., sizeof...(Is), counters[DIMC]);
				// logInfo("%td/%td", counters[DIMC], m_l);
			// if (R == REP::X)
				// logInfo("%td %g", counters[DIMC], std::norm(psi[counters[DIMC]]));

			if constexpr (std::is_same_v<Op, Identity>) res += std::norm(psi[counters[DIMC]]);
			else res += static_cast<Hamiltonian*>(this)->template call < Op >(InducedGrid::template pos<R>(counters[Is] + offset<R>(Is))...) * std::norm(psi[counters[DIMC]]);

		} while (!(...&& ((counters[reorder<Is>]++, counters[reorder<Is>] < (REP::X == R ? shape_l[reorder<Is>] : shape_t[reorder<Is>])) ? (counters[DIMC]++, false) : (counters[reorder<Is>] = 0, true))));

		// if (std::is_same_v<Op, KineticEnergy>)
			// _logMPI("%d %g %td %g", MPI::pID, res, counters[DIMC], test);

		return res * InducedGrid::template vol<R>();
	}


	template <REP R, class Op>
	inline double average()
	{
		return average_<R, Op>(indices);
	}
#pragma endregion Computations
	template <class Op>
	inline double getValue()
	{
		return 1.0;
	}
#pragma region Operations
	template <MODE M, REP R, class Op>
	inline auto operation()
	{
		if constexpr (std::is_same_v<Normalize, Op>)
		{
			auto res = average<R, Identity>();
			MPI::reduceImmediataly(&res);
			InducedGrid::multiplyArray(psi, sqrt(1.0 / res));
		}
		if constexpr (std::is_same_v<Symmetrize, Op>)
		{

		}
		if constexpr (std::is_same_v<AntiSymmetrize, Op>)
		{

		}
		if constexpr (std::is_same_v<Orthogonalize, Op>)
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
		// 		for (i = 0; i < m_l; i++)
		// 			wf[i] = wf[i] - amplits[lower] * ei[i];
		// 	}
		// }
		// auto res = average<R, Identity>(indices);
		// MPI::reduceImmediataly(&res);
		// InducedGrid::multiplyArray(psi, sqrt(1.0 / res));
		}
		// return res;
	}
#pragma endregion

#pragma region InitialConditions
	void setConstantValue(cxd val)
	{
		for (ind i = 0;i < m_l; i++) psi[i] = val;
	}

	template <class F, REP R, bool coords, uind ... Is>
	void add(F&& f, seq<Is...>)
	{
		ind dym;
		// _logMPI("my n0_o %td", n0_o);
		ind counters[DIMC + 1] = { 0 };
		do {
			if constexpr (coords)
				psi[counters[DIMC]] += f(InducedGrid::template pos<R>(counters[Is] + offset<R>(Is))...);
			else psi[counters[DIMC]] += f((counters[Is] + offset<R>(Is))...);
			// logInfo("%s %td %td", REP::X == R ? "X" : "P", counters[Is]...);

		} while (!(...&& ((counters[reorder<Is>]++, counters[reorder<Is>] < (REP::X == R ? shape_l[reorder<Is>] : shape_t[reorder<Is>]))
						  ? (counters[DIMC]++, false) : (counters[reorder<Is>] = 0, true))));
		// } while ((... ||
		// 		  !((counters[Is]++, counters[Is] < shape_l[Is]) ? (counters[DIMC]++, false) : (counters[Is] = 0, true))));
		// } while (!(... && ((counters[rev - Is]++, counters[rev - Is] < shape_l[rev - Is])
		// 				   ? (counters[DIMC]++, false)
		// 				   : (counters[rev - Is] = 0, true))));
						   // } while (!(
							   // ((counters[Is]++, counters[Is] < shape_l[Is]) ? (counters[DIMC]++, false) : (counters[Is] = 0, true)) && ... && true));
	}

	template <class F, REP R = REP::X>
	void addUsingCoordinateFunction(F&& f)
	{
		add<F, R, true>(std::forward<F>(f), indices);
	}
	template <class F, REP R = REP::X>
	void addUsingNodeFunction(F&& f)
	{
		add<F, R, false>(std::forward<F>(f), indices);
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