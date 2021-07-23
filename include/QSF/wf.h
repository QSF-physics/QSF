
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
	template <MODE M, REP R, class BO, class COMP>
	inline void compute(BO& bo, COMP&& c)
	{
		bo.template store <M, COMP>(1.0);
	}



	// ind constexpr mfix(ind Is)
	// {
	// 	return (Is == 0 ? local_start : 0);
	// }
	/* Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to
	switch x<->y for DIM>1 if working on opossite rep */
	template <ind Is>
	ind static constexpr swap01 = Is == 0 ? 1 : (Is == 1 ? 0 : Is);

	template <REP R, ind Is>
	ind static constexpr ofix = DIMC - 1 - ((R == REP::X) ? Is : swap01<Is>);
	/* When using mpi the first coordinate needs to be shifted
	NOTE: This is so when passing to functions accepting x,y,z...,
	internally the condition is different Is==0 for REP::X and Is==1 for REP::P*/

	ind constexpr mfix(ind Is) { return Is == 0 ? local_start : 0; };
	// {
	// 	if constexpr (R == REP::X) return DIMC - 1 - Is;
	// 	else return DIMC - 1 - (Is == 0 ? 1 : (Is == 1 ? 0 : Is));
	// }
	// static constexpr auto transposed_indices = swap_seq<n_seq_t<DIMC>>{};
	// //flip_seq_t<DIMC - 1, n_seq_t<DIMC>>;
	static constexpr auto indices = n_seq<DIMC>;
	// //flip_seq_t<DIMC - 1, n_seq_t<DIMC>>;
	// static constexpr ind rev = DIMC - 1;
	// // static constexpr n_seq_t<DIMC> indexes = n_seq<DIMC>;
	// template <REP R>
	// constexpr auto indices()
	// {	//FIXME: | !MPI::region
	// 	if constexpr (std::is_same_v<MPIStrategy, MPI::Single>)
	// 	{
	// 		if constexpr (R == REP::P && DIM > 1) return transposed_indices;
	// 		else return regular_indices;
	// 	}
	// 	if constexpr (!(std::is_same_v<MPIStrategy, MPI::Single>))
	// 	{
	// 		if ((R == REP::P && DIM > 1) && !MPI::region) return transposed_indices;
	// 		else return regular_indices;
	// 	}
	// 	// else
	// 	// {

	// 	// }
	// }

	template <MODE M, REP R, uind ... Is>
	inline void evolve_(double delta, seq<Is...>)
	{
		ind count[DIMC + 1] = { 0 };
		do {
			psi[count[DIMC]] *= static_cast<Hamiltonian*>(this)->template expOp < M, R>
				(delta * static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(count[Is] + mfix(Is))...));



			// if (R == REP::P)
			// {
			// 	psi[count[DIMC]] = { MPI::rID,
			// 	   static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(count[Is] + mfix(Is))...) };
			// }

			// if (R == REP::P)
				// _logMPI("%2d %2td [%2td %2td] (%2td %2td) [%2td %2td] %10g %10g %10g", MPI::rID, MPI::pID * local_m + count[DIMC],
						// (count[Is] + mfix(Is))..., (count[Is])...,
						// shape[Is]..., std::norm(psi[count[DIMC]]),
						// static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(count[Is] + mfix(Is))...), static_cast<Hamiltonian*>(this)->template expOp < M, R>
						// (delta * static_cast<Hamiltonian*>(this)->template operator() < R > (InducedGrid::template pos<R>(count[Is] + mfix(Is))...)));

			// static_cast<Hamiltonian*>(this)->
				// template expOp < M, R>(
					// delta * static_cast<Hamiltonian*>(this)->template call<std::conditional_t<R == REP::P, KineticEnergy, PotentialEnergy>>(InducedGrid::template pos<R>(count[Is] + mfix(Is))...));
		} while (!(...&& ((count[ofix<R, Is>]++, count[ofix<R, Is>] < shape[ofix<R, Is>])
						  ? (count[DIMC]++, false) : (count[ofix<R, Is>] = 0, true))));
	}
	template <MODE M, REP R>
	inline void evolve(double delta)
	{
		evolve_<M, R>(delta, indices);
	}


	template <REP R, class Op, uind ... Is>
	double average_(seq<Is...>)
	{
		// logInfo("%td %td %td %td", shape[0], shape[1], shape[2], sizeof...(Is));
		ind count[DIMC + 1]{ 0 };
		double res = 0.0;
		do {
			// logInfo("This will get added %s", typeid(Op).name());
			// if (R == REP::P)
				// logInfo("%2td %2td [%2td %2td] [%2td %2td] %td %td", Is..., count[Is]..., shape[Is]..., sizeof...(Is), count[DIMC]);
				// logInfo("%td/%td", count[DIMC], local_m);
			// if (R == REP::X)
				// logInfo("%td %g", count[DIMC], std::norm(psi[count[DIMC]]));

			if constexpr (std::is_same_v<Op, Identity>) res += std::norm(psi[count[DIMC]]);
			else res += static_cast<Hamiltonian*>(this)->template call < Op >(InducedGrid::template pos<R>(count[Is] + mfix(Is))...) * std::norm(psi[count[DIMC]]);

		} while (!(...&& ((count[ofix<R, Is>]++, count[ofix<R, Is>] < shape[ofix<R, Is>])
						  ? (count[DIMC]++, false) : (count[ofix<R, Is>] = 0, true))));

		// if (std::is_same_v<Op, KineticEnergy>)
			// _logMPI("%d %g %td %g", MPI::pID, res, count[DIMC], test);

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
		// 		for (i = 0; i < local_m; i++)
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
		for (ind i = 0;i < local_m; i++) psi[i] = val;
	}

	template <class F, REP R, bool coords, uind ... Is>
	void add(F&& f, seq<Is...>)
	{
		ind dym;
		// _logMPI("my local_start %td", local_start);
		ind count[DIMC + 1] = { 0 };
		do {
			if constexpr (coords)
				psi[count[DIMC]] += f(InducedGrid::template pos<R>(count[Is] + mfix(Is))...);
			else psi[count[DIMC]] += f((count[Is] + mfix(Is))...);
			// logInfo("%s %td %td", REP::X == R ? "X" : "P", count[Is]...);

		} while (!(...&& ((count[ofix<R, Is>]++, count[ofix<R, Is>] < shape[ofix<R, Is>])
						  ? (count[DIMC]++, false) : (count[ofix<R, Is>] = 0, true))));
		// } while ((... ||
		// 		  !((count[Is]++, count[Is] < shape[Is]) ? (count[DIMC]++, false) : (count[Is] = 0, true))));
		// } while (!(... && ((count[rev - Is]++, count[rev - Is] < shape[rev - Is])
		// 				   ? (count[DIMC]++, false)
		// 				   : (count[rev - Is] = 0, true))));
						   // } while (!(
							   // ((count[Is]++, count[Is] < shape[Is]) ? (count[DIMC]++, false) : (count[Is] = 0, true)) && ... && true));
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