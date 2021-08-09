#pragma region AvailableOperations
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
#pragma endregion AvailableOperations


template <class Hamiltonian, class BaseGrid, uind Components, class MPI_GC = typename BaseGrid::MPIGridComm, class MPIDiv = typename MPI_GC::MPIDivision, bool many = MPI_GC::many>
struct LocalGrid;

template <class MPIDivision>
struct GridPos;
template <>
struct GridPos<MPI::Slices>
{
	ind first;	//position of first local node in full array
	ind last;	//position of last local node in full array
};

template <class Hamiltonian, class BaseGrid, uind Components, class MPI_GC>
struct LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, false> : BaseGrid
{
	using BaseGrid::DIM;
	using BaseGrid::n;
	// using BaseGrid::inv_nn;
	using BaseGrid::inv_n;
	using BaseGrid::inv_2dx;
	using BaseGrid::nn;
	using BaseGrid::n2;
	using BaseGrid::xmin;
	using BaseGrid::dx, BaseGrid::dp;
	using MPIGridComm = MPI_GC;

	static constexpr auto indices = n_seq<DIM>;
	static constexpr bool hasAbsorber = std::is_base_of_v<AbsorberType, BaseGrid>;
	/* Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to
	switch x<->y for DIM>1 if working on opossite rep */
	template <uind Is>
	uind static constexpr swap01 = Is == 0 ? 1 : (Is == 1 ? 0 : Is);
	/* The N-dim loop passes indices in normal order, but increments in reverse */
	template <ind Is> ind static constexpr rev = DIM - 1 - Is;

	MPIGridComm mcomm; 			// MPI Grid Communicator
	const bool canLeaveTransposed;

	fftw_plan mpi_plans[2];


	cxd* psi_total = nullptr; 	// total ψ on m grid (used for saving)
	cxd* psi = nullptr;      	// local ψ(x,t) on m_l grid

	// Used in MPI calculations of currents
	cxd* row_before;
	cxd* row_after;
	borVec<DIM>* curr;

	const ind m = BaseGrid::m * Components;
	const double inv_m = BaseGrid::inv_m / Components;
	ind strides[DIM]; //The stride is the separation of consecutive elements along this dimension.
	ind m_l; //m local - total number of nodes per process

	GridPos<MPI::Slices> pos_lx;
	ind strides_lx[DIM];
	ind n_lx[DIM]; //MPI local n[]: local_n0,n1,n2,...

	// Different from shape_l[] if n0!=n1 and if FFTW_MPI_TRANSPOSED_OUT is used:
	GridPos<MPI::Slices> pos_lp;
	ind n_lp[DIM]; //P-space local MPI n[]: local_n1,n0,n2,...
	ind strides_lp[DIM];

	// const int transpose_f[2] = { 0, 0 };
	const unsigned transpose_f[2] = { FFTW_MPI_TRANSPOSED_OUT, FFTW_MPI_TRANSPOSED_IN };
	const int for_back[2] = { FFTW_FORWARD, FFTW_BACKWARD };

	LocalGrid(Section& settings) : BaseGrid(settings), mcomm(),
		canLeaveTransposed(DIM > 1 && mcomm.isMain),
		// mpiFFTW(mcomm.bounded[0]),
		m_l(m / MPI::rSize)
	{
		init();
		// mcomm.init(psi, m_l, n_lx, strides_lx);
		logSETUP("m %td m_l %td", m, m_l);
	}

	LocalGrid(BaseGrid g) : BaseGrid(g), mcomm(),
		canLeaveTransposed(DIM > 1 && mcomm.isMain),
		// mpiFFTW(mcomm.bounded[0]),
		m_l(m / MPI::rSize)
	{
		init();
		// mcomm.init(psi, m_l, n_lx, strides_lx);
		logSETUP("m %td m_l %td", m, m_l);
	}

	void gather()
	{
		if (!MPI::rID && psi_total == nullptr)
		{
			logALLOC("Allocating memory for psi_total (%td nodes)", m);
			psi_total = new cxd[m];
		}
		logMPI("Gathering " psi_symbol "... to address %p", psi_total);
		MPI_Gather(psi, m_l, MPI_CXX_DOUBLE_COMPLEX, psi_total, m_l, MPI_CXX_DOUBLE_COMPLEX, 0, MPI::rComm);
	}

	template <REP R, ind Is> ind constexpr abs_index(ind index) const noexcept
	{
		if constexpr (Is == 0)
		{
			if constexpr (REP::X == R) return index + pos_lx.first;
			else return index + pos_lp.first;
		}
		else return index;
	}
	template <REP R, ind dir> ind constexpr abs_centered_index(ind index) const noexcept
	{
		return abs_index<R, dir>(index) - n2[dir];
	}
	template <REP R, ind dir> ind constexpr neg_dist_from_edge(ind index) const noexcept
	{
		return abs_index<R, dir>(index) >= n2[dir] ? abs_index<R, dir>(index) - n[dir] + 1 : -abs_index<R, dir>(index);
	}


	template <REP R, uind dir> double abs_pos(ind index) noexcept
	{
		return BaseGrid::template pos<R, dir>(abs_index<R, dir>(index));
	}

	double abs2(ind index) { return std::norm(psi[index]); }

	template <REP R, uind dir> inline auto shape() const noexcept
	{
		return R == REP::X ? n_lx[dir] : n_lp[dir];
	}

	template <REP R, uind dir, uind dir_avg> inline auto bulk_shape() const noexcept
	{
		if constexpr (dir_avg != dir) return shape<R, dir>();
		else return shape<R, dir>() + ((dir_avg != 0 || MPI::rID == MPI::rSize - 1) ? -1 : 0);
	}

	template <uind dir, uind dir_avg> inline auto bulk_start() const noexcept
	{
		if constexpr (dir_avg != dir) return 0;
		else return ((dir_avg != 0 || MPI::rID == 0) ? 1 : 0);
	}

	template <REP R, ind dir> inline auto stride() const noexcept
	{
		return R == REP::X ? strides_lx[dir] : strides_lp[dir];
	}
	template <REP R> inline auto rowSize() const noexcept
	{
		return stride<R, 0>();
	}

	void initPositions()
	{
		pos_lx.first = MPI::rID * n_lx[0];
		pos_lx.last = (MPI::rID + 1) * n_lx[0] - 1;
		pos_lp.first = MPI::rID * n_lp[0];
		pos_lp.last = (MPI::rID + 1) * n_lp[0] - 1;
		// n0_e = n0_o + n0_l - 1;
	}
	template <REP R, uind ... Is>
	void initStrides(seq<Is...>)
	{
		ind* n_l = R == REP::X ? n_lx : n_lp;
		ind* strides_l = R == REP::X ? strides_lx : strides_lp;
		strides_l[DIM - 1] = 1;
		//Strides
		((strides_l[Is] = n_l[Is + 1] * strides_l[Is + 1]), ...);
	}
	template <REP R, uind ... Is>
	void initSizes(seq<Is...>)
	{
		ind* n_l = R == REP::X ? n_lx : n_lp;
		((n_l[Is] = n[Is]), ...);
		if constexpr (DIM > 1 && R == REP::P) if (canLeaveTransposed)
		{
			n_l[0] = n[1];
			n_l[1] = n[0];
		}
		n_l[0] = n[0] / MPI::rSize;
		// ((printf("n_l%s[%td]=%td\n", (R == REP::X ? "x" : "p"), Is, n_l[Is])), ...);
	}
	void initFFTW()
	{
		if (mcomm.isMain)
		{

			//make forward and backward mpi plans
			for (int i = 0; i < 2; i++)
				mpi_plans[i] =
				fftw_mpi_plan_many_dft(DIM, n, 1, 	//rank, dims,howmany
									   FFTW_MPI_DEFAULT_BLOCK, //block 
									   FFTW_MPI_DEFAULT_BLOCK, //tblock
									   reinterpret_cast<fftw_complex*>(psi),
									   reinterpret_cast<fftw_complex*>(psi),
									   MPI::rComm, for_back[i], MPI::plan_rigor | transpose_f[i]);
			reset();
		}
		// logSETUP("FFTW will be performed in one step");
		// logSETUP("Main Forward/backward plans for " psi_symbol " are initialized");
		// fftw_print_plan(transf_x2p);
		// fftw_print_plan(transf_p2x);
	}
	void init()
	{
		initSizes<REP::X>(n_seq<DIM>);
		initSizes<REP::P>(n_seq<DIM>);
		initStrides<REP::X>(rev_seq<DIM - 1>);
		initStrides<REP::P>(rev_seq<DIM - 1>);
		initPositions();
		testSizes();

		fftw_mpi_init();
		testFFTW();
		//Ready to alloc
		psi = (cxd*)fftw_malloc(sizeof(cxd) * m_l);


		/* Main region transforms all directions using MPI FFTW, others only use it to transform X.
		Transposing the first two dimensions in not all-dim cases would lead to problems */
		// ind nd[1] = { n_lx };
		initFFTW();
	}


#pragma region AccessOperators
	cxd& operator[](ind index) noexcept
	{
		return psi[index];
	}
	const cxd& operator[](ind index) const noexcept
	{
		return psi[index];
	}

	template <REP R, uind dim = 0, class I, class... Args>
	constexpr inline auto data_offset(I i, Args... args) // noexcept
	{
		if constexpr (DIM - 1 != dim) return i * stride<R, dim>() + data_offset<R, dim + 1>(args...);
		else  return i * stride<R, dim>();
	}

	template <REP R, class... I> cxd& operator()(I... i)
	{
		return psi[data_offset<R>(i...)];
	}
	template <REP R, class... I> const cxd& operator()(I... i) const
	{
		return psi[data_offset<R>(i...)];
	}

#pragma endregion AccessOperators

	//NOTE: Designed to work in X rep, FFTW_MPI_TRANSPOSED_OUT not taken into account
	template <REP R>
	void add(const cxd* psi_source, ind inp_n, double weight)
	{
		if (shape<R, 0>() == inp_n)
			for (ind i = 0; i < m_l; i++) psi[i] += weight * psi_source[i];
		else
		{
			ind m_down = (n[0] / 2 - inp_n / 2);
			ind m_up = m_down + inp_n - 1;
			for (ind i = m_down; i <= m_up; i++)
			{
				if (abs_pos<R>().first <= i && i <= abs_pos<R>().last)
				{
					if constexpr (DIM == 1)
						psi[i - abs_pos<R>().first] += weight * psi_source[i - m_down];
					else
					{
						for (ind j = m_down; j <= m_up; j++)
						{
							if constexpr (DIM == 2)
								psi[(i - abs_pos<R>().first) * n + j] += weight * psi_source[(i - m_down) * inp_n + (j - m_down)];
							else
							{
								for (ind k = m_down; k <= m_up; k++)
								{
									psi[(i - abs_pos<R>().first) * nn + j * n + k] += weight * psi_source[((i - m_down) * inp_n + (j - m_down)) * inp_n + (k - m_down)];
								}
							}
						}
					}
				}
			}
		}
	}

	ind localSize() { return m_l; }

	inline void reset() { resetArray(psi); }


	inline void normalizeAfterTwoFFT()
	{
		multiplyArray(psi, inv_m);
	}

	template <REP R>
	inline void fourier()
	{
	// Timings::measure::start("FFTW");
		static_assert(R == REP::X || R == REP::P, "Can only transform to X or P, unambigously.");
		constexpr DIMS back = R == REP::P ? 0 : 1;
		fftw_execute(mpi_plans[back]);
		if constexpr (back) normalizeAfterTwoFFT();
		// Timings::measure::stop("FFTW");
	}
	template <REP R, OPTIMS O>
	void precalc(double time)
	{
		static_cast<Hamiltonian*>(this)->_coupling.precalc(time);
	}
	void postCompute()
	{

	}

	//HACK: change second false to true after done with Dmitry!
	std::string save(std::string_view common_name = "",
					 DUMP_FORMAT df = { DIM, REP::X, true, true, true, true, false })
	{
		gather();
		if (!MPI::rID)
		{
			// logDUMPS("Dumping " psi_symbol " in %s rep", (R == REP::X ? "X" : "P"));
			// double dxORdp[DIM];
			// for (int i = 0;i < DIM; i++)
				// dxORdp[i] = mcomm.bounded[i] ? dx[i] : dp[i] / n[i];

			FILE* file = openPsi< AFTER<>, REP::X, DIM, IO_ATTR::WRITE>(common_name, 0, 0, true);

			writePsiBinaryHeader(file, n, dx, mcomm.bounded, df);
			fwrite(psi_total, sizeof(cxd), m, file);
			fclose(file);
		}
		return std::string(file_path);
	}

#pragma region Computations
	//TODO: if no match here pass to derived class
	template <MODE M, REP R, class BO, class COMP>
	inline void compute(BO& bo, COMP&& c)
	{
		bo.template store <M, COMP>(1.0);
	}


	/* When using mpi the first coordinate needs to be shifted
	NOTE: This is so when passing to functions accepting x,y,z...*/
	template <MODE M, REP R, uind ... Is>
	inline void evolve_(double delta, seq<Is...>)
	{
		ind counters[DIM + 1] = { 0 };
		do {
			// logInfo("pos [%g %g] val: %g delta: %g rep: %d", abs_pos<R, Is>(counters[Is])..., static_cast<Hamiltonian*>(this)->template operator() < R, Is... > (abs_pos<R, Is>(counters[Is])...), delta, ind(R));
			// if (mcomm.isMain)
			psi[counters[DIM]] *=
				static_cast<Hamiltonian*>(this)->template expOp < M, R>(delta * static_cast<Hamiltonian*>(this)->template operator() < R, Is... > (abs_pos<R, Is>(counters[Is])...));

			if constexpr (hasAbsorber && R == REP::X && !MPIGridComm::many)
			{
				psi[counters[DIM]] *= BaseGrid::template absorb<Is...>(delta, neg_dist_from_edge<R, Is>(counters[Is])...);
				// logWarning("Calling with %g", BaseGrid::template operator() < Is... > (abs_centered_index<R, Is>(counters[Is])...));
			}

		} while (!(...&&
				   ((counters[rev<Is>]++, counters[rev<Is>] < shape<R, rev<Is>>())
					? (counters[DIM]++, false)
					: (counters[rev<Is>] = 0, true))
				   ));
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
		ind counters[DIM + 1]{ 0 };
		double res = 0.0;
		do {
			res += abs2(counters[DIM]) * static_cast<Hamiltonian*>(this)->template
				call < Op, Is... >(abs_pos<R, Is>(counters[Is])...);

		} while (!(...&&
				   ((counters[rev<Is>]++, counters[rev<Is>] < shape<R, rev<Is>>())
					? (counters[DIM]++, false)
					: (counters[rev<Is>] = 0, true))
				   ));

		return res * BaseGrid::template vol<R>();
	}
	template <REP R, class Op>
	inline double average()
	{
		return average_<R, Op>(indices);
	}

	template <REP R, uind dir_avg, class Op, uind ... Is>
	double average_der_(seq<Is...>)
	{
		// logInfo("%td %td %td %td", shape_l[0], shape_l[1], shape_l[2], sizeof...(Is));
		ind counters[DIM + 1]{ 0 };
		counters[dir_avg] = 0;
		double res = 0.0;
		do {
			res += abs2(counters[DIM]) *
				(static_cast<Hamiltonian*>(this)->template
				 call < Op, Is... >(abs_pos<R, Is>(counters[Is] + (Is == dir_avg ? 1 : 0))...)
				 - static_cast<Hamiltonian*>(this)->template
				 call < Op, Is... >(abs_pos<R, Is>(counters[Is] + (Is == dir_avg ? -1 : 0))...));

		} while (!(...&&
				   ((counters[rev<Is>]++, counters[rev<Is>] < shape<R, rev<Is>>())
					? (counters[DIM]++, false)
					: (counters[rev<Is>] = 0, true))
				   ));

		return res * BaseGrid::template vol<R>() * inv_2dx[dir_avg];
	}
	template <REP R, uind dir, class Op>
	inline double average_der()
	{
		return average_der_<R, dir, Op>(indices);
	}


	template <REP R, uind dir>
	inline void update_curr(ind counter) //imag(conj(psi))*∂psi
	{
		if (counter - stride<R, dir>() < 0 || counter + stride<R, dir>() > m_l)
			curr[counter][dir] = 0.0;
		else
			curr[counter][dir] = inv_2dx[dir]
			* imag(conj(psi[counter]) * (psi[counter + stride<R, dir>()] - psi[counter - stride<R, dir>()]));

		// logWarning("%g", curr[counter][dir]);
	}
	template <REP R, uind ... Is>
	inline void current_map_(seq<Is...>)
	{
		ind counter = 0;
		while (counter < m_l)
		{
			(update_curr<R, Is>(counter), ...);
			counter++;
		}

		if (MPI::rSize > 1) //edge cases for MPI
		{
			if (MPI::pID < MPI::rSize - 1) // Recieve from higher MPI::rID
			{
				for (ind j = 0; j < rowSize<R>(); j++)
					curr[j][0] = inv_2dx[0]
					* imag(conj(psi[j]) * (row_after[j] - psi[j - rowSize<R>()]));
			}
			if (MPI::pID > 0) // Recieved from lower MPI::rID
			{
				for (ind j = 0; j < rowSize<R>(); j++)
					curr[j][0] = inv_2dx[0]
					* imag(conj(psi[j]) * (psi[j + rowSize<R>()] - row_before[j]));
			}
		}
	}
	template <REP R>
	inline void current_map()
	{
		if (curr == nullptr) //TODO: move to prepare
		{
			curr = (borVec<DIM>*)malloc(sizeof(borVec<DIM>) * m_l);
			if (MPI::rSize > 1)
			{
				row_after = new cxd[rowSize<R>()];
				row_before = new cxd[rowSize<R>()];
			}
		}
		getNeighbouringNodes<R>();
		current_map_<R>(indices);
	}
	template <REP R>
	void getNeighbouringNodes()
	{
		if (MPI::rSize > 1)
		{
			int tag = 1;
			if (MPI::rID > 0) // Send up
				MPI_Send(psi, int(rowSize<R>()), MPI_CXX_DOUBLE_COMPLEX, MPI::rID - 1, tag, MPI::rComm);
			if (MPI::rID < MPI::rSize - 1)	// Recieve
				MPI_Recv(row_after, int(rowSize<R>()), MPI_CXX_DOUBLE_COMPLEX, MPI::rID + 1, tag, MPI::rComm, &MPI::status);
			tag = 2;
			if (MPI::rID < MPI::rSize - 1)	// Send down
				MPI_Send(&(psi[(shape<R, 0>() - 1) * rowSize<R>()]), int(rowSize<R>()), MPI_CXX_DOUBLE_COMPLEX, MPI::rID + 1, tag, MPI::rComm);
			if (MPI::rID > 0)			// Recieve 
				MPI_Recv(row_before, int(rowSize<R>()), MPI_CXX_DOUBLE_COMPLEX, MPI::rID - 1, tag, MPI::rComm, &MPI::status);
		}
	}
	template <REP R, class FLUX_TYPE, uind ... Is>
	inline double flux_(seq<Is...>)
	{
		ind counters[DIM + 1]{ 0 };
		double res = 0.0;
		do {
			// logWarning("%g", FLUX_TYPE::border((counters[Is] - n2[Is])...)[0]);
			// logWarning("%g %g %g", curr[counters[DIM]][Is]..., FLUX_TYPE::border((counters[Is] - n2[Is])...)[0], res);
			res += curr[counters[DIM]]
				* FLUX_TYPE::border((abs_centered_index<R, Is>(counters[Is]))...);

		} while (!(...&&
				   ((counters[rev<Is>]++, counters[rev<Is>] < shape<R, rev<Is>>())
					? (counters[DIM]++, false)
					: (counters[rev<Is>] = 0, true))
				   ));
		// return res * BaseGrid::template vol<R>();
		return res * BaseGrid::template vol<R>() / dx[0];
	}
	template <REP R, class FLUX_TYPE, uind ... Is>
	inline double flux()
	{
		return flux_<R, FLUX_TYPE>(indices);
	}

	template <class Op>
	inline double getValue()
	{
		return static_cast<Hamiltonian*>(this)->_coupling.template getValue<Op>();
	}
#pragma endregion Computations

#pragma region Operations
	template <REP R, uind ... Is>
	inline void maskRegion_(seq<Is...>)
	{
		ind counters[DIM + 1] = { 0 };
		do
		{
			psi[counters[DIM]] *= BaseGrid::template mask<Is...>(
				(mcomm.bounded[Is] ? neg_dist_from_edge<R, Is>(counters[Is]) : 0)...);

		} while (!(...&&
				   ((counters[rev<Is>]++, counters[rev<Is>] < shape<R, rev<Is>>())
					? (counters[DIM]++, false)
					: (counters[rev<Is>] = 0, true))
				   ));
	}
	inline void maskRegion()
	{
		maskRegion_<REP::X>(indices);
	}

	template <MODE M, REP R, class Op>
	inline auto operation()
	{
		if constexpr (std::is_same_v<Normalize, Op>)
		{
			auto res = average<R, Identity>();
			MPI::reduceImmediataly(&res);
			multiplyArray(psi, sqrt(1.0 / res));
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
		// multiplyArray(this->psi, sqrt(1.0 / res));
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
		ind counters[DIM + 1] = { 0 };
		do {
			if constexpr (coords) psi[counters[DIM]] += f(abs_pos<R, Is>(counters[Is])...);
			else psi[counters[DIM]] += f((abs_index<R, Is>(counters[Is]))...);

		} while (!(...&&
				   ((counters[rev<Is>]++, counters[rev<Is>] < shape<R, rev<Is>>())
					? (counters[DIM]++, false)
					: (counters[rev<Is>] = 0, true))
				   ));
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

	void load(std::string input_path)
	{
		if (!MPI::pID)
		{
			// auto out = PsiFile(M_FROM, )
			FILE* fin = fopen_with_check<IO_ATTR::READ>(input_path.c_str());
			//openPsi<AFTER<>, REP::X, DIM, IO::ATTR::READ>(name, 0, 0, true);
			bool is_complex = readPsiBinaryHeader<DIM>(fin);

			if (is_complex) fread(psi, sizeof(cxd), m, fin);
			else for (ind i = 0; i < m; i++) fread(&psi[i], sizeof(double), 1, fin);
			closeFile(fin);

			double qnorm = 0.0;
			for (ind i = 0; i < m; i++) qnorm += std::norm(psi[i]);
			logSETUP("State " psi_symbol "_%d loaded with norm %g", 0, BaseGrid::template vol<REP::X>() * qnorm);
		}
	}


#pragma region ArrayOperations
	//BUNCH OF HELPFUL FUNCTIONS
	inline void resetArray(cxd* array)
	{
		for (ind i = 0; i < m_l; i++) array[i] = 0.0;
	}
	inline void multiplyArray(cxd* array, const double mult)
	{
		// logInfo("Multiplying array by %g", mult);
		for (ind i = 0; i < m_l; i++) array[i] *= mult;
	}
	inline void multiplyArray(cxd* array, const double* mult)
	{
		for (ind i = 0; i < m_l; i++)
			array[i] *= mult[i];
	}
	inline void copyArray(cxd* to, cxd* from)
	{
		for (ind i = 0; i < m_l; i++)
			to[i] = from[i];
	}

	inline void copyFullArray(cxd* to, cxd* from)
	{
		for (ind i = 0; i < m; i++)
			to[i] = from[i];
	}

	inline void premultiplyFullArray(cxd* array, const double mult)
	{
		for (ind i = 0; i < m; i++)
			array[i] *= mult;
	}

	inline cxd scalarProduct(cxd* to, cxd* from)
	{
		cxd over = 0.;
		for (ind i = 0; i < m_l; i++)
			over += conj(to[i]) * from[i];
		return over;
	}

	inline cxd scalarProduct(double* to, cxd* from)
	{
		cxd over = 0.;
		for (ind i = 0; i < m_l; i++)
			over += to[i] * from[i];
		return over;
	}
	inline cxd overlap(cxd* psi_to, cxd* psi_from, double mult)
	{
		cxd over = scalarProduct(psi_to, psi_from);
		return over * mult;
	}
	inline cxd overlap(double* psi_to, cxd* psi_from, double mult)
	{
		cxd over = scalarProduct(psi_to, psi_from);
		return over * mult;
	}

#pragma endregion ArrayOperations

#pragma endregion InitialConditions

#pragma region Tests

	void testSizes()
	{
		if (!(n[0] % MPI::rSize == 0)) logWarning("Grid length n[0] (%td) should be divisible by the number of MPI region processes (%d) (FFTW reg)", n[0], MPI::rSize);

		if constexpr (DIM == 1)
		{
			if (!(m % (MPI::rSize * MPI::rSize) == 0))
				logWarning("Grid size m=n[0] (%td) should be divisible by the number of MPI region processes squared (%d) (FFTW req)", n[0], MPI::pSize * MPI::pSize);
		}
		else
		{
			if (!((m / n[0]) % (MPI::rSize) == 0))
				logWarning("Row size (m/n[0]=%td) should be divisible by the number of MPI processes (%d) (FFTW req)", m / n[0], MPI::rSize);
		}
		logSETUP("Setting up " psi_symbol "(x,t) arrays and plans...");
		for (DIMS i = 0; i < DIM; i++) logSETUP("Sizes n_lx[%d]=%td", i, n_lx[i]);
		for (DIMS i = 0; i < DIM; i++) logSETUP("Sizes n_lp[%d]=%td", i, n_lp[i]);
		for (DIMS i = 0; i < DIM; i++) logSETUP("Sizes strides_lx[%d]=%td", i, strides_lx[i]);
		for (DIMS i = 0; i < DIM; i++) logSETUP("Sizes strides_lp[%d]=%td", i, strides_lp[i]);
	}
	void testFFTW()
	{
		if constexpr (DIM == 1)
		{
			/* We need to deal with 1D case seperately, see:
			http://fftw.org/doc/Basic-and-advanced-distribution-interfaces.html#Basic-and-advanced-distribution-interfaces */
			ind ret_n0_lx_i, ret_pos_lx_first_i;
			ind ret_n0_lx_o, ret_pos_lx_first_o;
			ind ret_m_l = fftw_mpi_local_size_1d(n[0], MPI::rComm,
												 FFTW_FORWARD, MPI::plan_rigor,
												 &ret_n0_lx_i, &ret_pos_lx_first_i,
												 &ret_n0_lx_o, &ret_pos_lx_first_o);

			if (m_l != ret_m_l) logWarning("Expected m_l [%td] different from the one return by fftw_mpi_local_size_1d [%td]! Abort.", m_l, ret_m_l);
			if (n_lx[0] != ret_n0_lx_i) logWarning("Expected n_lx[0] [%td] different from the one return by fftw_mpi_local_size_1d [%td] (input)! Abort.", n_lx[0], ret_n0_lx_i);
			if (n_lx[0] != ret_n0_lx_o) logWarning("Expected n_lx[0] [%td] different from the one return by fftw_mpi_local_size_1d [%td] (output)! Abort.", n_lx[0], ret_n0_lx_o);
			//TODO: do other tests
		}
		else
		{
			ind ret_n0_lx_i, ret_pos_lx_first_i;
			ind ret_m_l = fftw_mpi_local_size(DIM, n, MPI::rComm,
											  &ret_n0_lx_i, &ret_pos_lx_first_i);
			if (m_l != ret_m_l) logWarning("Expected m_l [%td] different from the one return by fftw_mpi_local_size_1d [%td]! Abort.", m_l, ret_m_l);
			//TODO: do other tests
		}
	}
#pragma endregion Tests
};



