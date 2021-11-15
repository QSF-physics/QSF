namespace QSF
{

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

	template <DIMS D>
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
		cxd* row_before = nullptr;
		cxd* row_after = nullptr;
		borVec<DIM>* curr = nullptr;

		const ind m = BaseGrid::m * Components;
		const double inv_m = BaseGrid::inv_m / Components;
		ind strides[DIM]; //The stride is the separation of consecutive elements along this dimension.
		ind m_l; //m local - total number of nodes per process

		GridPos<MPI::Slices> pos_lx;
		ind n_lx[DIM]; //MPI local n[]: local_n0,n1,n2,...
		ind strides_lx[DIM];

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
			logSETUP("m %td m_l %td", m, m_l);
		}

		LocalGrid(BaseGrid g) : BaseGrid(g), mcomm(),
			canLeaveTransposed(DIM > 1 && mcomm.isMain),
			// mpiFFTW(mcomm.bounded[0]),
			m_l(m / MPI::rSize)
		{
			init();
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

		void scatter()
		{
			logMPI("Scattering " psi_symbol "... to address %p", psi);
			MPI_Scatter(psi_total, m_l, MPI_CXX_DOUBLE_COMPLEX, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, 0, MPI::rComm);
		}

		template <REP R, ind Is> ind constexpr abs_index(ind index) const noexcept
		{
			if constexpr (Is == 0)
			{
				if constexpr (R == REP::X) return index + pos_lx.first;
				else return index + pos_lp.first;
			}
			else return index;
		}
		template <REP R, ind dir> ind constexpr abs_centered_index(ind index) const noexcept
		{
			return abs_index<R, dir>(index) - n2[dir];
		}
		//only for X rep
		template <ind dir> ind constexpr neg_dist_from_edge(ind index) const noexcept
		{
			return abs_index<REP::X, dir>(index) >= n2[dir]
				? abs_index<REP::X, dir>(index) - n[dir] + 1
				: -abs_index<REP::X, dir>(index);
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
		}
		template <uind ... Is>
		void initAbsStrides(seq<Is...>)
		{
			strides[DIM - 1] = 1;
			((strides[Is] = n[Is + 1] * strides[Is + 1]), ...);
		}
		template <REP R, uind ... Is>
		void initLocalStrides(seq<Is...>)
		{
			ind* n_l = R == REP::X ? n_lx : n_lp;
			ind* strides_l = R == REP::X ? strides_lx : strides_lp;
			strides_l[DIM - 1] = 1;
			((strides_l[Is] = n_l[Is + 1] * strides_l[Is + 1]), ...);
		}
		template <REP R, uind ... Is>
		void initLocalSizes(seq<Is...>)
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
			initLocalSizes<REP::X>(n_seq<DIM>);
			initLocalSizes<REP::P>(n_seq<DIM>);
			initAbsStrides(rev_seq<DIM - 1>);
			initLocalStrides<REP::X>(rev_seq<DIM - 1>);
			initLocalStrides<REP::P>(rev_seq<DIM - 1>);
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
		template <uind dim = 0, class I, class... Args>
		constexpr inline auto abs_data_offset(I i, Args... args) // noexcept
		{
			if constexpr (DIM - 1 != dim) return i * strides[dim] + abs_data_offset<dim + 1>(args...);
			else  return i * strides[dim];
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
			// logInfo("%g", inv_m);
			multiplyArray(psi, inv_m);
		}

		template <REP R>
		inline void fourier()
		{
			// logInfo("FT %d", int(R));
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
				// if (((counters[rev<Is>] == 511) && ...))
					// logInfo("pos [0,0] = [%g %g] val: %g delta: %g rep: %d", abs_pos<R, Is>(counters[Is])..., static_cast<Hamiltonian*>(this)->template operator() < R, Is... > (abs_pos<R, Is>(counters[Is])...), delta, ind(R));
					// if (mcomm.isMain)
				psi[counters[DIM]] *=
					static_cast<Hamiltonian*>(this)->template expOp <M, R>(delta * static_cast<Hamiltonian*>(this)->template operator() < R, Is... > (abs_pos<R, Is>(counters[Is])...));

				if constexpr (hasAbsorber && R == REP::X && !MPIGridComm::many)
				{
					psi[counters[DIM]] *= BaseGrid::template absorb<Is...>(delta, neg_dist_from_edge<Is>(counters[Is])...);
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
			// logInfo("evolve %d", int(R));
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
			// logInfo("norm^2: %g", res * BaseGrid::template vol<R>());
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


		template <uind dir>
		inline void update_curr(ind counter) //imag(conj(psi))*∂psi
		{
			cxd tmp;
			if ((counter - stride<REP::X, dir>() > 0) && (counter + stride<REP::X, dir>() < m_l)) //Boundries
			{
				// fprintf(stderr, "%d %td %td %p %p %p\n", MPI::pID, stride<R, dir>(), counter, &psi[counter], &psi[counter + stride<R, dir>()], &psi[counter - stride<R, dir>()]);
				tmp = psi[counter + stride<REP::X, dir>()] - psi[counter - stride<REP::X, dir>()];
				curr[counter][dir] = inv_2dx[dir] * imag(conj(psi[counter]) * tmp);
			}
			else curr[counter][dir] = 0.0;
		}
		template <uind ... Is>
		inline void current_map_(seq<Is...>)
		{
			ind counter = 0;
			// printf("SIZES: [%td %td]\n", m, m_l);
			while (counter < m_l)
			{
				(update_curr<Is>(counter), ...);
				counter++;
			}

			if (MPI::rSize > 1) //edge cases for MPI
			{
				if (MPI::rID < MPI::rSize - 1) // Recieve from higher MPI::rID
				{
					ind lastRow = m_l - rowSize<REP::X>();
					for (ind j = 0; j < rowSize<REP::X>(); j++)
					{
						// fprintf(stderr, "%d %td %td %p %p %p\n", MPI::pID, rowSize<R>(), j, &psi[j], &psi[j - rowSize<R>()], &row_after[j]);
						curr[lastRow + j][0] = inv_2dx[0]
							* imag(conj(psi[lastRow + j]) * (row_after[j] - psi[lastRow - rowSize<REP::X>() + j]));
					}
				}
				if (MPI::rID > 0) // Recieved from lower MPI::rID
				{
					for (ind j = 0; j < rowSize<REP::X>(); j++)
					{
						// fprintf(stderr, "%d %td %td %p %p %p\n", MPI::pID, rowSize<R>(), j, &psi[j], &psi[j + rowSize<R>()], &row_before[j]);
						curr[j][0] = inv_2dx[0]
							* imag(conj(psi[j]) * (psi[j + rowSize<REP::X>()] - row_before[j]));
					}
				}
			}
		}

		inline borVec<DIM>* current_map()
		{
			//TODO: make it fail in P rep
			//TODO: move to prepare 
			if (curr == nullptr) curr = new borVec<DIM>[m_l];

			getNeighbouringNodes(psi, rowSize<REP::X>(), shape<REP::X, 0>());
			current_map_(indices);
			return curr;
		}

		void getNeighbouringNodes(cxd* source, ind rowSize, ind rowCount)
		{

			//TODO: change to MPI_window
			if (MPI::rSize > 1)
			{
				if (row_after == nullptr) row_after = new cxd[rowSize];
				if (row_before == nullptr) row_before = new cxd[rowSize];

				if (MPI::rID > 0) // Send up
					MPI_Send(source, rowSize, MPI_CXX_DOUBLE_COMPLEX, MPI::rID - 1, 13, MPI::rComm);
				if (MPI::rID < MPI::rSize - 1)	// Recieve
					MPI_Recv(row_after, rowSize, MPI_CXX_DOUBLE_COMPLEX, MPI::rID + 1, 13, MPI::rComm, &MPI::status);

				if (MPI::rID < MPI::rSize - 1)	// Send down
					MPI_Send(&(source[(rowCount - 1) * rowSize]), rowSize, MPI_CXX_DOUBLE_COMPLEX, MPI::rID + 1, 7, MPI::rComm);
				if (MPI::rID > 0)			// Recieve 
					MPI_Recv(row_before, rowSize, MPI_CXX_DOUBLE_COMPLEX, MPI::rID - 1, 7, MPI::rComm, &MPI::status);
			}
		}

		template <REP R, class FLUX_TYPE, uind ... Is>
		inline double flux_(seq<Is...>)
		{
			ind counters[DIM + 1]{ 0 };
			double res = 0.0;
			do {
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
		template <REP R, class FLUX_TYPE>
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
		template <uind ... Is>
		inline void maskRegion_(seq<Is...>)
		{
			ind counters[DIM + 1] = { 0 };
			do
			{
				psi[counters[DIM]] *= BaseGrid::template mask<Is...>(
					(mcomm.bounded[Is] ? neg_dist_from_edge<Is>(counters[Is]) : -n2[Is])...);

			} while (!(...&&
					   ((counters[rev<Is>]++,
						 counters[rev<Is>] < shape<REP::X, rev<Is>>())
						? (counters[DIM]++, false)
						: (counters[rev<Is>] = 0, true))
					   ));
		}
		// REP::X only!
		inline void maskRegion()
		{
			maskRegion_(indices);
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
				//TODO: raise not implemented!
			}
			if constexpr (std::is_same_v<AntiSymmetrize<3_D>, Op>)
			{
				if constexpr (DIM < 3) return; //TODO:: log error
				else
				{
					if (MPI::rSize > 1)
					{
						if (psi_total == nullptr)
							psi_total = new cxd[m];
					   // MPI_Gather(psi, m_l, MPI_CXX_DOUBLE_COMPLEX, psi_total, m_l, MPI_CXX_DOUBLE_COMPLEX, 0, MPI::rComm);
						MPI_Allgather(psi, m_l, MPI_CXX_DOUBLE_COMPLEX, psi_total, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI::rComm);
						ind regular;
						for (ind i = 0; i < n_lx[0]; i++)
						{
							for (ind j = 0; j < n_lx[1]; j++)
							{
								for (ind k = 0; k < n_lx[2]; k++)
								{
									regular = data_offset<REP::X>(i, j, k);
									// psi[regular] = psi_total[abs_data_offset(i + pos_lx.first, j, k)]; //Already there
									psi[regular] += psi_total[abs_data_offset(j, k, i + pos_lx.first)];
									psi[regular] += psi_total[abs_data_offset(k, i + pos_lx.first, j)];
									psi[regular] -= psi_total[abs_data_offset(j, i + pos_lx.first, k)];
									psi[regular] -= psi_total[abs_data_offset(i + pos_lx.first, k, j)];
									psi[regular] -= psi_total[abs_data_offset(k, j, i + pos_lx.first)];
								}
							}
						}
					}
					else
					{
						ind regular;
						ind regular2;
						ind regular3;
						ind switched;
						ind switched2;
						ind switched3;
						for (ind i = 0; i < n_lx[0]; i++)
						{
							for (ind j = 0; j < i; j++)
							{
								for (ind k = 0; k < j; k++)
								{
									regular = data_offset<REP::X>(i, j, k);
									regular2 = data_offset<REP::X>(j, k, i);
									regular3 = data_offset<REP::X>(k, i, j);
									switched = data_offset<REP::X>(j, i, k);
									switched2 = data_offset<REP::X>(i, k, j);
									switched3 = data_offset<REP::X>(k, j, i);
									psi[regular] += psi[regular2] + psi[regular3]
										- psi[switched] - psi[switched2] - psi[switched3];
									psi[regular2] = psi[regular];
									psi[regular3] = psi[regular];
									psi[switched] = -psi[regular];
									psi[switched2] = -psi[regular];
									psi[switched3] = -psi[regular];
								}
							}
							for (ind j = 0; j < n_lx[1]; j++)
							{
								psi[data_offset<REP::X>(i, j, j)] = 0.0;
								psi[data_offset<REP::X>(j, i, j)] = 0.0;
								psi[data_offset<REP::X>(j, j, i)] = 0.0;
							}
						}
					}
				}
			}
			if constexpr (std::is_same_v<AntiSymmetrize<2_D>, Op>)
			{
				if constexpr (DIM < 2) return;
				else
				{
					if constexpr (DIM == 3_D) //Partial anti-symmetrization as in Neon
					{
						ind regular;
						ind switched;
						for (ind i = 0; i < n_lx[0]; i++)
						{
							for (ind j = 0; j < n_lx[1]; j++)
							{
								for (ind k = 0; k < j; k++)
								{
									regular = data_offset<REP::X>(i, j, k);
									switched = data_offset<REP::X>(i, k, j);
									psi[regular] -= psi[switched];
									psi[switched] = -psi[regular];
								}
								psi[data_offset<REP::X>(i, j, j)] = 0.0;
							}
						}
					}
					else
					{
						if (psi_total == nullptr) psi_total = new cxd[m];
					   // MPI_Gather(psi, m_l, MPI_CXX_DOUBLE_COMPLEX, psi_total, m_l, MPI_CXX_DOUBLE_COMPLEX, 0, MPI::rComm);
						MPI_Allgather(psi, m_l, MPI_CXX_DOUBLE_COMPLEX, psi_total, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI::rComm);
						ind regular;
						for (ind i = 0; i < n_lx[0]; i++)
						{
							for (ind j = 0; j < n_lx[1]; j++)
							{
								regular = data_offset<REP::X>(i, j);// (i * n + j) * n + k;
								psi[regular] -= psi_total[abs_data_offset(j, i + pos_lx.first)];
							}
						}
					}
				}
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
	#pragma endregion InitialConditions

	#pragma region IO
		void _load(IO::path input_path, std::string ext = IO::psi_ext)
		{
			input_path += ext;

			MPI_File fh;
			DUMP_FORMAT df;
			MPI_File_open(MPI::rComm, input_path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
			logIO("[_load] Opening [%s] file %s", "rb", input_path.c_str());

				//TODO: verify
			MPI_File_read(fh, &df.dim, 1, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
			if (df.dim != DIM)
			{
				logError("Dimensionality of the input file is not the same as the one choosen for evolution");
			}
			MPI_File_read(fh, &df.rep, 1, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);

			MPI_File_read(fh, &df.unnormalized, 1, MPI_CXX_BOOL, MPI_STATUS_IGNORE);
			MPI_File_read(fh, &df.initial_wf_subtracted, 1, MPI_CXX_BOOL, MPI_STATUS_IGNORE);
			int size;
			MPI_File_read(fh, &size, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
			MPI_File_read(fh, n, df.dim, MPI_INT64_T, MPI_STATUS_IGNORE);
			MPI_File_read(fh, dx, df.dim, MPI_DOUBLE, MPI_STATUS_IGNORE);
			// MPI_File_read(fh, mcomm.bounded, df.dim, MPI_CXX_BOOL, MPI_STATUS_IGNORE);
			// IO::writePsiBinaryHeader(file, n, dx, mcomm.bounded, df);

			int offset = sizeof(DIMS) + sizeof(REP) + sizeof(bool) + sizeof(bool);
			offset += sizeof(int) + DIM * (sizeof(ind) + sizeof(double) + sizeof(bool));
			MPI_File_set_view(fh, offset + MPI::rID * m_l * sizeof(cxd), MPI_CXX_DOUBLE_COMPLEX,
							  MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
			// MPI_File_read_all(fh, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
			// MPI_File_seek(fh, offset + MPI::rID * m_l, MPI_SEEK_SET);
			MPI_File_read(fh, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
			MPI_File_close(&fh);
			// if (MPI::region == 0)
			{
				auto res = average<REP::X, Identity>();
				MPI::reduceImmediataly(&res);
				__logMPI("pID %d, State " psi_symbol "_%d loaded with norm %g\n", MPI::pID, 0, res);
			}
			for (ind i = 0; i < m_l; i++)
			{
				psi[i] = { std::abs(psi[i]),0 };
			}
		// if (MPI::region == region_index)
		// {
		// 	if (!MPI::rID) //only node can load the file
		// 	{
		// 		if (psi_total == nullptr)
		// 		{
		// 			logALLOC("Allocating memory for psi_total (%td nodes) loading file %s", m, input_path.c_str());
		// 			psi_total = new cxd[m];
		// 		}
		// 			// auto out = PsiFile(M_FROM, )
		// 		FILE* fin = fopen(input_path.c_str(), "rb");
		// 		//IO::fopen_with_check<IO::IO_ATTR::READ>(input_path.c_str());
		// 		//openPsi<AFTER<>, REP::X, DIM, IO::ATTR::READ>(name, 0, 0, true);
		// 		bool is_complex = IO::readPsiBinaryHeader<DIM>(fin);

		// 		if (is_complex) fread(psi_total, sizeof(cxd), m, fin);
		// 		else for (ind i = 0; i < m; i++) fread(&psi_total[i], sizeof(double), 1, fin);
		// 		fclose(fin);


		// 		double qnorm = 0.0;
		// 		for (ind i = 0; i < m; i++) qnorm += std::norm(psi_total[i]);
		// 		logSETUP("State " psi_symbol "_%d loaded with norm %g", 0, BaseGrid::template vol<REP::X>() * qnorm);
		// 	}
		// 	scatter();
		// }
		}

		void croossOut()
		{
			ind ip4 = n[0] / 2 + (ind)(50.0 / dx[0]);
			ind ip3 = n[0] / 2 - (ind)(50.0 / dx[0]);
			logInfo("crossOut");
			ind intemp2 = ip4 - ip3;
			ind intemp;
			// y-direction
			for (ind i = ip3;i < ip4;i++) {
				double temp = 0.5 * (1. + cos(2. * pi * (i - ip3) / intemp2));
				for (ind j = 0;j < n[0];j++)
				{
					if (i >= pos_lx.first && i <= pos_lx.last)
					{
						intemp = (i - pos_lx.first) * n[0] + j;
						psi[intemp] = temp * psi[intemp];
						// logInfo("x-crossout %td %td %td %g", intemp, i, j, psi[intemp]);
					}
				}
			}
			// x-direction
			for (ind j = ip3;j < ip4;j++) {
				double temp = 0.5 * (1. + cos(2. * pi * (j - ip3) / intemp2));
				for (ind i = 0;i < n[0];i++) {
					if (i >= pos_lx.first && i <= pos_lx.last)
					{
						intemp = (i - pos_lx.first) * n[0] + j;
						psi[intemp] = temp * psi[intemp];
						// logInfo("y-crossout %td %td %td %g", intemp, i, j, psi[intemp]);
					}
				}
			}
		}

		void orthogonalizeWith(IO::path input_path, std::string ext = IO::psi_ext)
		{
			input_path += ext;

			MPI_File fh;
			MPI_File_open(MPI::rComm, input_path.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
			logIO("[orthogonalizeWith] Opening [%s] file %s", "rb", input_path.c_str());


			int offset = sizeof(DIMS) + sizeof(REP) + sizeof(bool) + sizeof(bool);
			offset += sizeof(int) + DIM * (sizeof(ind) + sizeof(double) + sizeof(bool));
			MPI_File_set_view(fh, offset + MPI::rID * m_l * sizeof(cxd), MPI_CXX_DOUBLE_COMPLEX,
							  MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
			// MPI_File_read_all(fh, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
			// MPI_File_seek(fh, offset + MPI::rID * m_l, MPI_SEEK_SET);
			cxd* tmp = new cxd[m_l];
			MPI_File_read(fh, tmp, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
			MPI_File_close(&fh);
			cxd scalarProduct = 0.0;
			for (ind i = 0; i < m_l; i++)
				scalarProduct += conj(tmp[i]) * psi[i];
			scalarProduct *= BaseGrid::template vol<REP::X>();
			MPI::reduceImmediataly(&scalarProduct);
			logIO("|<psi_0|psi>|^2 =%g", std::norm(scalarProduct));
			double _norm = 0;
			for (ind i = 0; i < m_l; i++)
			{
				psi[i] -= scalarProduct * tmp[i];
				_norm += std::norm(psi[i]);
			}
			_norm *= BaseGrid::template vol<REP::X>();
			MPI::reduceImmediataly(&_norm);
			logIO("After orthogonalization |<psi|psi>|^2 =%g", _norm);
			delete[] tmp;
		}

		void load(IO::path input_path, std::string ext = IO::psi_ext)
		{
			_load(input_path, ext);
		}
		void restore() {
			logWarning("Restored the wf");
			_load("backup");
		}

		IO::path _save(IO::path input_path,
					   DUMP_FORMAT df,
					   std::string ext = IO::psi_ext)
		{
			input_path += ext;

			MPI_File fh;
			// MPI_File_delete(input_path.c_str(), MPI_INFO_NULL);
			MPI_File_open(MPI::rComm, input_path.c_str(),
						  MPI_MODE_CREATE | MPI_MODE_WRONLY,
						  MPI_INFO_NULL, &fh);

			logIO("Opening [%s] file %s", "wb", input_path.c_str());

			if (MPI::rID == 0)
			{
				//TODO: save step at which it is for restart
				MPI_File_write(fh, &df.dim, 1, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
				MPI_File_write(fh, &df.rep, 1, MPI_UNSIGNED_CHAR, MPI_STATUS_IGNORE);
				MPI_File_write(fh, &df.unnormalized, 1, MPI_CXX_BOOL, MPI_STATUS_IGNORE);
				MPI_File_write(fh, &df.initial_wf_subtracted, 1, MPI_CXX_BOOL, MPI_STATUS_IGNORE);
				int size = true ? sizeof(cxd) : sizeof(double);
				MPI_File_write(fh, &size, 1, MPI_INT32_T, MPI_STATUS_IGNORE);
				MPI_File_write(fh, n, df.dim, MPI_INT64_T, MPI_STATUS_IGNORE);
				MPI_File_write(fh, dx, df.dim, MPI_DOUBLE, MPI_STATUS_IGNORE);
				MPI_File_write(fh, mcomm.bounded, df.dim, MPI_CXX_BOOL, MPI_STATUS_IGNORE);
				// IO::writePsiBinaryHeader(file, n, dx, mcomm.bounded, df);
			}
			int offset = sizeof(df.dim) + sizeof(df.rep) + sizeof(df.rep) + sizeof(df.initial_wf_subtracted);
			offset += sizeof(int) + df.dim * (sizeof(ind) + sizeof(double) + sizeof(bool));
			MPI_File_set_view(fh, offset + MPI::rID * m_l * sizeof(cxd), MPI_CXX_DOUBLE_COMPLEX,
							  MPI_CXX_DOUBLE_COMPLEX, "native", MPI_INFO_NULL);
			MPI_File_write(fh, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
			// MPI_File_seek(fh, offset + MPI::rID * m_l, MPI_SEEK_SET);
			// MPI_File_write(fh, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_STATUS_IGNORE);
			MPI_File_close(&fh);

		// else
		// {
		// 	gather();
		// 	if ((!mainRegionOnly && !MPI::rID) || (mainRegionOnly && !MPI::pID))
		// 	{
		// 		FILE* file;
		// 		logWarning("%s =?= %s %d", path.parent_path().c_str(), IO::target_path.c_str(), path.parent_path().compare(IO::target_path));

		// 		if ((path.has_parent_path() && path.parent_path().compare(IO::target_path)) || path.has_extension()) file = fopen(path.c_str(), "wb");
		// 		else
		// 			file = IO::openPsi< AFTER<>, REP::X, DIM, IO::IO_ATTR::WRITE>(path.c_str(), 0, 0, true);
		// 	   // logDUMPS("Dumping " psi_symbol " in %s rep", (R == REP::X ? "X" : "P"));

		// 		IO::writePsiBinaryHeader(file, n, dx, mcomm.bounded, df);
		// 		fwrite(psi_total, sizeof(cxd), m, file);
		// 		fclose(file);
		// }

			return std::filesystem::current_path() / input_path;
		}

		template <uind ... Is>
		void phaseShiftTrick(seq<Is...>)
		{
			ind counters[DIM + 1]{ 0 };
			bool odd = 0;
			do {
				odd = (counters[Is] + ...) % 2;
				if (odd)
					psi[counters[DIM]] *= -1;
			} while (!(...&&
					   ((counters[rev<Is>]++, counters[rev<Is>] < shape<REP::X, rev<Is>>())
						? (counters[DIM]++, false)
						: (counters[rev<Is>] = 0, true))
					   ));
		}
		IO::path save(IO::path path = "", DUMP_FORMAT df = { .dim = DIM })
		{
			if (df.rep == REP::P)
			{
				logIO("Applying phaseShiftTrick to move to centered P representation");
				phaseShiftTrick(indices);
				logIO("Transforming to P representation");
				fourier<REP::P>();
			}
			auto ret = _save(path, df);
			if (df.rep == REP::P)
			{
				logIO("Transforming back to X representation");
				fourier<REP::X>();
				phaseShiftTrick(indices);
			}
			return ret;
		}
		void _backup(const ind step, std::string ext = IO::psi_ext)
		{
			if (!MPI::pID) {
				FILE* f = IO::fopen_with_check("backup.info", "wb");
				fwrite(&step, sizeof(ind), 1, f);
				fclose(f);
			}
			_save("backup", { .dim = DIM }, ext);
			logInfo("Backup made at step %td", step);
		}
		void backup(const ind step)
		{
			_backup(step);
		}
		void remove_backups(std::string ext = IO::psi_ext)
		{
			if (!MPI::rID) //main node of each region removes
				std::filesystem::remove("backup" + ext);
		}
		void snapshot(std::string extra_info, DUMP_FORMAT df = { .dim = DIM })
		{
			static ind x_count = 0;
			static ind p_count = 0;
			if (df.rep == REP::X)
			{
				save("snapshot_" + std::to_string(x_count) + extra_info, df);
				x_count++;
			}
			else
			{
				save("snapshot_" + std::to_string(p_count) + extra_info, df);
				p_count++;
			}
		}
	#pragma endregion IO

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



}