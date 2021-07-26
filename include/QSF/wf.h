


// template <class Hamiltonian, class GridBase, size_t Components>
// struct WF : Grid <GridBase, Components>
// {
// 	// using H = Hamiltonian; //Derived type!
// 	using G = Grid <GridBase, Components>;
// 	using G::DIMC;
// 	using G::DIM;
// 	using G::kin_scale;
// 	using G::psi;
// 	using G::abs2;

// 	WF(Section& settings) : G(settings) {
// 		logInfo("WF init");
// 	}

// 	WF(GridBase g) : G(g) {
// 		logInfo("WF init");
// 	}

// #pragma region Computations
// 	//TODO: if no match here pass to derived class
// 	template <MODE M, REP R, class BO, class COMP>
// 	inline void compute(BO& bo, COMP&& c)
// 	{
// 		bo.template store <M, COMP>(1.0);
// 	}
// 	/* The N-dim loop passes indices in normal order, but increments in reverse */
// 	template <ind Is>
// 	ind static constexpr rev = DIMC - 1 - Is;

// 	/* When using mpi the first coordinate needs to be shifted
// 	NOTE: This is so when passing to functions accepting x,y,z...,
// 	internally the condition is different Is==0 for REP::X and Is==1 for REP::P*/


// 	template <MODE M, REP R, uind ... Is>
// 	inline void evolve_(double delta, seq<Is...>)
// 	{
// 		ind counters[DIMC + 1] = { 0 };
// 		do {
// 			G::psi[counters[DIMC]] *= static_cast<Hamiltonian*>(this)->template expOp < M, R>(delta * static_cast<Hamiltonian*>(this)->template operator() < R > (G::template abs_pos<R, Is>(counters[Is])...));

// 			// if (R == REP::P)
// 			// {
// 			// 	psi[counters[DIMC]] = { MPI::rID,
// 			// 	   static_cast<Hamiltonian*>(this)->template operator() < R > (G::template pos<R>(counters[Is] + offset<R>(Is))...) };
// 			// }

// 			// if (R == REP::P)
// 			// 	_logMPI("%2d %2td [%2td %2td] (%2td %2td) [%2td %2td] %10g %10g %10g", MPI::rID, MPI::pID * m_l + counters[DIMC],
// 			// 			(counters[Is] + offset<R>(Is))..., (counters[Is])...,
// 			// 			shape_l[Is]..., std::norm(psi[counters[DIMC]]),
// 			// 			static_cast<Hamiltonian*>(this)->template operator() < R > (G::template pos<R>(counters[Is] + offset<R>(Is))...), static_cast<Hamiltonian*>(this)->template expOp < M, R>
// 			// 			(delta * static_cast<Hamiltonian*>(this)->template operator() < R > (G::template pos<R>(counters[Is] + offset<R>(Is))...)));

// 			// static_cast<Hamiltonian*>(this)->
// 				// template expOp < M, R>(
// 					// delta * static_cast<Hamiltonian*>(this)->template call<std::conditional_t<R == REP::P, KineticEnergy, PotentialEnergy>>(G::template pos<R>(counters[Is] + offset<R>(Is))...));
// 		} while (
// 			!(...&&
// 			  ((counters[rev<Is>]++, counters[rev<Is>] < G::template shape<R>(rev<Is>)
// 				? (counters[DIMC]++, false)
// 				: (counters[rev<Is>] = 0, true)))));
// 	}
// 	template <MODE M, REP R>
// 	inline void evolve(double delta)
// 	{
// 		evolve_<M, R>(delta, G::indices);
// 	}

// 	template <REP R, class Op, uind ... Is>
// 	double average_(seq<Is...>)
// 	{
// 		// logInfo("%td %td %td %td", shape_l[0], shape_l[1], shape_l[2], sizeof...(Is));
// 		ind counters[DIMC + 1]{ 0 };
// 		double res = 0.0;
// 		do {
// 			// logInfo("This will get added %s", typeid(Op).name());
// 			// if (R == REP::P)
// 				// logInfo("%2td %2td [%2td %2td] [%2td %2td] %td %td", Is..., counters[Is]..., shape_l[Is]..., sizeof...(Is), counters[DIMC]);
// 				// logInfo("%td/%td", counters[DIMC], m_l);
// 			// if (R == REP::X)
// 				// logInfo("%td %g", counters[DIMC], std::norm(psi[counters[DIMC]]));

// 			if constexpr (std::is_same_v<Op, Identity>)
// 				res += G::abs2(counters[DIMC]);
// 			else res += static_cast<Hamiltonian*>(this)->template call < Op >(G::template abs_pos<R, Is>(counters[Is])...) * G::abs2(counters[DIMC]);

// 		} while (!(...&& ((counters[rev<Is>]++, counters[rev<Is>] < G::template shape<R>(rev<Is>)) ? (counters[DIMC]++, false) : (counters[rev<Is>] = 0, true))));

// 		// if (std::is_same_v<Op, KineticEnergy>)
// 			// _logMPI("%d %g %td %g", MPI::pID, res, counters[DIMC], test);

// 		return res * G::template vol<R>();
// 	}


// 	template <REP R, class Op>
// 	inline double average()
// 	{
// 		return average_<R, Op>(G::indices);
// 	}
// #pragma endregion Computations
// 	template <class Op>
// 	inline double getValue()
// 	{
// 		return 1.0;
// 	}
// #pragma region Operations
// 	template <MODE M, REP R, class Op>
// 	inline auto operation()
// 	{
// 		if constexpr (std::is_same_v<Normalize, Op>)
// 		{
// 			auto res = average<R, Identity>();
// 			MPI::reduceImmediataly(&res);
// 			G::multiplyArray(G::psi, sqrt(1.0 / res));
// 		}
// 		if constexpr (std::is_same_v<Symmetrize, Op>)
// 		{

// 		}
// 		if constexpr (std::is_same_v<AntiSymmetrize, Op>)
// 		{

// 		}
// 		if constexpr (std::is_same_v<Orthogonalize, Op>)
// 		{

// 		// if (state > 0)
// 		// {
// 		// 	int lower = 0;
// 		// 	([&] {
// 		// 		if (lower < state)
// 		// 			amplits[lower] = _CO_PROJ<IdentityOperator<REP::X | REP::P>>::template calc<R, opt, integral_constant<size_t, Args>>();
// 		// 		lower++;
// 		// 	 }(), ...);

// 		// 	MPI::reduceImmediataly(amplits, amplits_size);
// 		// 	for (lower = 0; lower < state; lower++)
// 		// 	{
// 		// 		cplxd* ei = wf->states[lower];
// 		// 		for (i = 0; i < m_l; i++)
// 		// 			wf[i] = wf[i] - amplits[lower] * ei[i];
// 		// 	}
// 		// }
// 		// auto res = average<R, Identity>(indices);
// 		// MPI::reduceImmediataly(&res);
// 		// G::multiplyArray(this->psi, sqrt(1.0 / res));
// 		}
// 		// return res;
// 	}
// #pragma endregion

// #pragma region InitialConditions
// 	void setConstantValue(cxd val)
// 	{
// 		for (ind i = 0;i < G::m_l; i++) G::psi[i] = val;
// 	}

// 	template <class F, REP R, bool coords, uind ... Is>
// 	void add(F&& f, seq<Is...>)
// 	{
// 		// _logMPI("my n0_o %td", n0_o);
// 		ind counters[DIMC + 1] = { 0 };
// 		do {
// 			if constexpr (coords)
// 				G::psi[counters[DIMC]] += f(G::template abs_pos<R, Is>(counters[Is])...);
// 			else G::psi[counters[DIMC]] += f((G::template abs_index<R, Is>(counters[Is]))...);
// 			// logInfo("%s %td %td", REP::X == R ? "X" : "P", counters[Is]...);

// 		} while (!(...&& ((counters[rev<Is>]++, counters[rev<Is>] < G::template shape<R>(rev<Is>)) ? (counters[DIMC]++, false) : (counters[rev<Is>] = 0, true))));
// 		// } while ((... ||
// 		// 		  !((counters[Is]++, counters[Is] < shape_l[Is]) ? (counters[DIMC]++, false) : (counters[Is] = 0, true))));
// 		// } while (!(... && ((counters[rev - Is]++, counters[rev - Is] < shape_l[rev - Is])
// 		// 				   ? (counters[DIMC]++, false)
// 		// 				   : (counters[rev - Is] = 0, true))));
// 						   // } while (!(
// 							   // ((counters[Is]++, counters[Is] < shape_l[Is]) ? (counters[DIMC]++, false) : (counters[Is] = 0, true)) && ... && true));
// 	}

// 	template <class F, REP R = REP::X>
// 	void addUsingCoordinateFunction(F&& f)
// 	{
// 		add<F, R, true>(std::forward<F>(f), this->indices);
// 	}
// 	template <class F, REP R = REP::X>
// 	void addUsingNodeFunction(F&& f)
// 	{
// 		add<F, R, false>(std::forward<F>(f), this->indices);
// 	}

// 	void addFromFile()
// 	{

// 	}

// #pragma endregion InitialConditions

// };

