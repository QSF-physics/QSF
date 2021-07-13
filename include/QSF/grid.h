template <class BaseGrid, size_t Components, class MPIStrategy = typename BaseGrid::MPIStrategy, class MPIGrids = typename BaseGrid::MPIGrids>
struct Grid;

template <class BaseGrid, size_t Components>
struct Grid<BaseGrid, Components, MPI::Slices> : BaseGrid
{
	// 	static_assert(std::is_base_of_v< CoordinateSystem, CS>,
// 				  "First GridType template argument must be type derived from CoordinateSystem");

// 	// static_assert(is_base_of_v<MPI::Grid, MR>,
// 				//   "Third GridType template argument must be type derived from MPI::Grid");

// 	static_assert((!std::is_same_v<MPI::Rods, MPIStrategy>) || (DIM >= DIMS::D3), "Rod MPI strategy not supported for dimensionality < 3");
	using BaseGrid::DIM;
	using MPISol = MPIGrid<typename BaseGrid::MPIGrids, BaseGrid::D, MPI::Slices>;
	MPISol sol;
	using BaseGrid::absorber;
	using BaseGrid::m;
	using BaseGrid::n;
	using BaseGrid::inv_m;
	using BaseGrid::nn;
	cxd* psi_total = nullptr; // total ψ on m grid (used for X and P)
	cxd* psi;      			// local ψ(x,t) on local_m grid
	cxd* psi_copy;      		// backup ψ(x,t) on local_m grid
	cxd* psi_acc;      		// accumulated ψ(x,t) on local_m grid
	ind local_n;
	ind local_m;
	ind local_start;
	ind local_end;

	ind ns[DIM]; //dimensions
	ind strides[DIM];
	ind local_ns[DIM]; //dimensions
	ind local_strides[DIM];

	fftw_plan mpi_plans[2];
	bool mpiFFTW;
	fftw_plan extra_plans[2 * DIM];
	int extra_plans_count;

	void gather()
	{
		if (!MPI::rID && psi_total == nullptr)
		{
			logALLOC("Allocating memory for psi_total");
			psi_total = new cxd[m];
		}
		logMPI("Gathering " psi_symbol "...");
		MPI_Gather(psi, local_m, MPI_CXX_DOUBLE_COMPLEX, psi_total, local_m, MPI_CXX_DOUBLE_COMPLEX, 0, MPI::rComm);
	}
	cxd& operator[](ind index)
	{
		return psi[index];
	}

	Grid() : BaseGrid()
	{
		fftw_mpi_init();
		local_n = n / MPI::rSize; //Expected split
		strides[DIM - 1] = 1;
		ns[DIM - 1] = n;

		for (int i = DIM - 2; i >= 0; i--)
		{
			ns[i] = n; //Size of each dimension
			strides[i] = ns[i + 1] * strides[i + 1];
		}

		logSETUP("Setting up " psi_symbol "(x,t) arrays and plans...");
		logSETUP("Grid size n (%td) will split into local_n (%td) by rSize (%d)",
				 n, local_n, MPI::rSize);

		if constexpr (DIM == 1)
		{
			/* We need to deal with 1D case seperately, see:
			http://fftw.org/doc/Basic-and-advanced-distribution-interfaces.html#Basic-and-advanced-distribution-interfaces */
			ind local_n2, local_start2; //Used in 1D only
			local_n2 = local_n;
			local_m = fftw_mpi_local_size_1d(n, MPI::rComm, FFTW_FORWARD, MPI::plan_rigor, &local_n, &local_start, &local_n2, &local_start2);
			// local_m = fftw_mpi_local_size_1d(n, MPI::rComm, FFTW_BACKWARD, plan_rigor, &local_n, &local_start, &local_n2, &local_start2);
		}
		else local_m = fftw_mpi_local_size(DIM, ns, MPI::rComm, &local_n, &local_start);
		psi = (cxd*)fftw_malloc(sizeof(cxd) * local_m);
		// We should use non-MPI FFTW routines when transform doesn't involve X (is a free coord)

		// Transposing the first two dimensions in not all-dim cases would lead to problems
		mpiFFTW = sol.bounded[0];
		const bool canLeaveTransposed = DIM > 1 && MPI::isMain;
		int transpose_f[2] = { FFTW_MPI_TRANSPOSED_OUT, FFTW_MPI_TRANSPOSED_IN };
		int for_back[2] = { FFTW_FORWARD, FFTW_BACKWARD };
		ind nd[1] = { n };

/* Main region transforms all directions using MPI FFTW, others only use it to transform X */
		if (mpiFFTW) for (int i = 0; i < 2; i++) //make forward and backward mpi plans
			mpi_plans[i] = fftw_mpi_plan_many_dft(MPI::isMain ? DIM : 1, //rank
												  MPI::isMain ? ns : nd, //dims
												  MPI::isMain ? 1 : m / ns[0], //howmany
												  FFTW_MPI_DEFAULT_BLOCK, //block 
												  FFTW_MPI_DEFAULT_BLOCK, //tblock
												  reinterpret_cast<fftw_complex*>(psi),
												  reinterpret_cast<fftw_complex*>(psi),
												  MPI::rComm, for_back[i], MPI::plan_rigor | (canLeaveTransposed ? transpose_f[i] : 0));

/* Non-mpi fftw (extra) plans for non-main region (region!=0)
transforms are coupled together if involve consecutive directions (case for 1-3D)*/
		ns[0] = local_n;
		int fftw_rank = 0;
		int start_dim = -1;
		extra_plans_count = 0;
		if (!MPI::isMain) for (int i = 1; i <= DIM; i++)
		{
			if (i == DIM || !sol.bounded[i])
			{
				if (fftw_rank)
				{
					fftw_iodim64* dims = new fftw_iodim64[fftw_rank];
					for (int j = 0; j < fftw_rank; j++)
					{
						dims[j] = { ns[start_dim + j], strides[start_dim + j],strides[start_dim + j] };
						__logMPI("region %d start_dim %d j %d ns[start_dim+j]=%td strides[start_dim+j]=%td\n", MPI::region, start_dim, j, ns[start_dim + j], strides[start_dim + j]);
					}

					//define a loop over all howmany_dims (free) dimensions
					int index = 0;
					int howmany_rank = DIM - fftw_rank;
					fftw_iodim64* howmany_dims = new fftw_iodim64[howmany_rank];
					for (int j = 0; j < DIM;j++)
					{
						if (j < start_dim || (j >= start_dim + fftw_rank && j < DIM))
						{
							howmany_dims[index] = { ns[j], strides[j], strides[j] };
							index++;
							__logMPI("region %d start_dim %d index %d j %d ns[j]=%td dist[j]=%td\n", MPI::region, start_dim, index, j, ns[j], strides[j]);
						}
					};
					for (int k = 0; k < 2; k++)
						extra_plans[extra_plans_count + DIM * k] = fftw_plan_guru64_dft(
							fftw_rank, dims,
							howmany_rank, howmany_dims,
							reinterpret_cast<fftw_complex*>(psi),
							reinterpret_cast<fftw_complex*>(psi),
							for_back[k], MPI::plan_rigor);

					delete[] dims;
					delete[] howmany_dims;
					extra_plans_count++;
					__logMPI("region %d mpiFFTW %d extra_plans_count %d dim %d/%d fftw_rank %d howmany_rank %d\n", MPI::region, mpiFFTW, extra_plans_count, start_dim, i - 1, fftw_rank, howmany_rank);
				}
				fftw_rank = 0;
			}
			else if (sol.bounded[i])
			{
				if (start_dim == -1)start_dim = i;
				fftw_rank++;
			}
		}

		__logMPI("region %d mpiFFTW %d extra_plans_count %d\n", MPI::region, mpiFFTW, extra_plans_count);


		local_end = local_start + local_n - 1;

		_logMPI("Region %d MPI::pID: %d got %td nodes or %td rows [start row: %td end row: %td]", MPI::region, MPI::pID, local_m, local_n, local_start, local_end);
		logSETUP("Psi and forward/backward plans initialized");
		// fftw_print_plan(transf_x2p);
		// fftw_print_plan(transf_p2x);
	}
		//Designed to work in X rep, FFTW_MPI_TRANSPOSED_OUT not taken into account
	void add(const cxd* psi_source, ind inp_n, double weight)
	{
		if (n == inp_n)
			for (ind i = 0; i < local_m; i++) psi[i] += weight * psi_source[i];
		else
		{
			ind m_down = (n / 2 - inp_n / 2);
			ind m_up = m_down + inp_n - 1;
			for (ind i = m_down; i <= m_up; i++)
			{
				if (local_start <= i && i <= local_end)
				{
					if constexpr (DIM == 1)
						psi[i - local_start] += weight * psi_source[i - m_down];
					else
					{
						for (ind j = m_down; j <= m_up; j++)
						{
							if constexpr (DIM == 2)
								psi[(i - local_start) * n + j] += weight * psi_source[(i - m_down) * inp_n + (j - m_down)];
							else
							{
								for (ind k = m_down; k <= m_up; k++)
								{
									psi[(i - local_start) * nn + j * n + k] += weight * psi_source[((i - m_down) * inp_n + (j - m_down)) * inp_n + (k - m_down)];
								}
							}
						}
					}
				}
			}
		}
	}
	inline void reset()
	{
		resetArray(psi);
	}

	void backup()
	{
		for (ind i = 0; i < local_m; i++)
		{
			psi_copy[i] = psi[i];
			psi_acc[i] = 0.0;
		}
	}

	void restore()
	{
		for (ind i = 0; i < local_m; i++) psi[i] = psi_copy[i];
	}


	void accumulate(double coeff)
	{
		for (ind i = 0; i < local_m; i++) psi_acc[i] += coeff * psi[i];
	}

	void collect(double coeff)
	{
		for (ind i = 0; i < local_m; i++)
			psi[i] = coeff * psi[i] + psi_acc[i];
	}
	inline void normalizeAfterTwoFFT()
	{
		switch (sol.freeCoord)
		{
		case FREE_COORD::NO:
			premultiplyArray(psi, inv_m); break;
		case FREE_COORD::X:
		case FREE_COORD::Y:
		case FREE_COORD::Z:
			premultiplyArray(psi, 1.0 / nn); break;
		case FREE_COORD::XY:
		case FREE_COORD::YZ:
		case FREE_COORD::XZ:
			premultiplyArray(psi, 1.0 / n); break;
		default: break;
		}
	}

	template <REP R>
	inline void fourier()
	{
		Timings::measure::start("FFTW");
		static_assert(R == REP::X || R == REP::P, "Can only transform to X or P, unambigously.");
		constexpr size_t shift = R == REP::P ? 0 : 1;
		if (mpiFFTW) fftw_execute(mpi_plans[shift]);
		for (int i = 0; i < extra_plans_count; i++)
			fftw_execute(extra_plans[i + DIM * shift]);
		if constexpr (shift) premultiplyArray(psi, inv_m * (pow(n, DIM - sol.boundedCoordDim)));
		Timings::measure::stop("FFTW");

	}
	//BUNCH OF HELPFUL FUNCTIONS
	inline void copyFullArray(cxd* to, cxd* from)
	{
		for (ind i = 0; i < m; i++)
			to[i] = from[i];
	}
	inline void copyArray(cxd* to, cxd* from)
	{
		for (ind i = 0; i < local_m; i++)
			to[i] = from[i];
	}
	inline void premultiplyFullArray(cxd* array, const double mult)
	{
		for (ind i = 0; i < m; i++)
			array[i] *= mult;
	}
	inline void premultiplyArray(cxd* array, const double mult)
	{
		for (ind i = 0; i < local_m; i++)
			array[i] *= mult;
	}
	inline void multiplyArrayBy(cxd* array, const double* mult)
	{
		for (ind i = 0; i < local_m; i++)
			array[i] *= mult[i];
	}
	inline void resetArray(cxd* array)
	{
		for (ind i = 0; i < local_m; i++)
			array[i] = 0.;
	}
	inline cxd scalarProduct(cxd* to, cxd* from)
	{
		cxd over = 0.;
		for (ind i = 0; i < local_m; i++)
			over += conj(to[i]) * from[i];
		return over;
	}

	inline cxd scalarProduct(double* to, cxd* from)
	{
		cxd over = 0.;
		for (ind i = 0; i < local_m; i++)
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
};




template <class BaseGrid, size_t Components>
struct Grid<BaseGrid, Components, MPI::Slices, MPI::Multi> : Grid<BaseGrid, Components, MPI::Slices>
{
	using base = Grid<BaseGrid, Components, MPI::Slices>;
	using base::DIM;
	using base::sol;
	using base::psi;
	using base::m;
	using base::n;
	using base::xmin;
	using base::dx;
	using base::nn;
	using base::local_n;
	using base::absorber;
	using base::local_start;
	using base::local_end;
	static inline constexpr cxd zero = { 0.0, 0.0 };
	ind slice_n;
	ind slice_start;
	ind slice_n2;
	ind slice_start2;
	ind slice_m;
	fftw_plan transf_slice[BaseGrid::DIM];

	MPI_Win lessFreeWin = nullptr;
	MPI_Win moreFreeWin = nullptr;

	cxd* slice;
	ind sliceSize;
	int boxesCount;
	int nCAP;
	bool CAPnodesPerP;

	// void initSlice();
	void initWindow()
	{
		if (sol.moreFree) MPI_Win_create(psi, m * sizeof(cxd), sizeof(cxd), MPI_INFO_NULL, sol.moreFree, &moreFreeWin);
		if (sol.lessFree) MPI_Win_create(NULL, 0, sizeof(cxd), MPI_INFO_NULL, sol.lessFree, &lessFreeWin);
	}
	Grid() :base()
	{
		// absorber.correct(n, L);
		// nCAP = absorber.nCAP;
		// boxesCount = DIM < 3 ? 1 : n / nCAP;
		// CAPnodesPerP = local_n / nCAP;
		initSlice();
		initWindow();
	}
	~Grid()
	{
		// destroySlice();
	}

	//See: http://www.fftw.org/fftw3_doc/Load-balancing.html#Load-balancing
	void test()
	{
		// logTestFatal(n % MPI::rSize == 0, "Grid length n (%td) should be divisible by the number of MPI region processes (%d) (FFTW reg)", n, MPI::rSize);
		if constexpr (BaseGrid::DIM == 1)
		{
			// logTestFatal(m % (MPI::rSize * MPI::rSize) == 0, "Grid length n (%td) should be divisible by the number of MPI region processes squared (%d) (FFTW req)", n, MPI::pSize * MPI::pSize);
		}
		else
		{
			// logTestFatal((m / n) % (MPI::rSize) == 0, "Row size (m/n=%td) should be divisible by the number of MPI processes (%d) (FFTW req)", m / n, MPI::rSize);
		}

		// logTestFatal(n >= Im::n && n >= Src::n, "Grid length n (%td) should be larger than IM grid length Im::n (%td) and previous mode grid length Src::n (%td)", n, Im::n, Src::n);
	}



	void initSlice()
	{
		/* TODO: Here we decide how the wf will be sliced.
		1) SIZE: In 1D slices are not cut, hence sliceSize == n.
		In other cases size of the slice should be constant for all regions sliceSize == nCAP*m/n.
		2) We always want to split the job among all the processes in the region, hence the slices
		will be cut along the X-direction, whenever we are transfering a direction != X.
		When transfering X direction we cut along Y direction as we need to perform fftw along X. */
		if (sol.boundedCoordDim == DIM) return; //Main region doesn't operate with slices

		sliceSize = DIM == 1 ? n : (DIM == 2 ? nn : absorber.nCAP * nn);
		ind howmany = DIM == 1 ? 1 : (DIM == 2 ? n : absorber.nCAP * n);
		ind dims[1] = { n };
		slice = (cxd*)fftw_malloc(sizeof(cxd) * sliceSize / MPI::rSize);
		transf_slice[0] = fftw_mpi_plan_many_dft(1, dims, howmany,
												 FFTW_MPI_DEFAULT_BLOCK, //block 
												 FFTW_MPI_DEFAULT_BLOCK, //tblock
												 reinterpret_cast<fftw_complex*>(slice),
												 reinterpret_cast<fftw_complex*>(slice),
												 MPI::rComm, FFTW_FORWARD, MPI::plan_rigor);
		if (transf_slice[0] == NULL) logInfo("fftw_mpi_plan... (slice) returned NULL");

		ind ns[DIM];
		ns[0] = local_n;
		ind strides[DIM];
		strides[DIM - 1] = 1;
		for (int i = 1; i < DIM; i++)
		{
			//HACK: very dirty
			if constexpr (DIM == 2) ns[1] = n;
			else if constexpr (DIM == 3)
			{
				ns[1] = (i == 1 ? n : nCAP);
				ns[2] = (i == 1 ? nCAP : n);
			}

			for (int j = DIM - 2; j >= 0; j--)
				strides[j] = ns[j + 1] * strides[j + 1];

			fftw_iodim64 dims[1] = { n , strides[i], strides[i] };
			fftw_iodim64 howmany_dims[DIM - 1];
			int index = 0;
			for (int j = 0; j < DIM; j++)
			{
				if (j != i)
				{
					howmany_dims[index] = { ns[j], strides[j],strides[j] };
					index++;
				}
			}
			printf("i=%d ns[i]=%td %td\n", i, ns[i], strides[i]);
			transf_slice[i] = fftw_plan_guru64_dft(
				1, dims, DIM - 1, howmany_dims, //rank, dims, howmany_rank, howmany_dims
				reinterpret_cast<fftw_complex*>(slice),
				reinterpret_cast<fftw_complex*>(slice),
				FFTW_FORWARD, MPI::plan_rigor);
		}
	}


	void destroyWindow()
	{
		if (sol.moreFree) MPI_Win_free(&moreFreeWin);
		if (sol.lessFree) MPI_Win_free(&lessFreeWin);
	}

	// (0,0,0) = (bottom, left, front)
	inline void getXBox(int rank, int boxIndex)
	{
		if constexpr (DIM == 1)
		{
			if (local_start < nCAP) //bottom
			{
				int size = Min(int(local_n), int(nCAP) - int(local_start));
				MPI_Get(&slice[0], size, MPI_CXX_DOUBLE_COMPLEX, rank,
						0, size, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
			if (local_end > n - 1 - nCAP) //top
			{
				int start = Max(0, int(n) - int(nCAP) - int(local_start));
				MPI_Get(&slice[start], local_n - start, MPI_CXX_DOUBLE_COMPLEX, rank,
						start, local_n - start, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
		}
		if constexpr (DIM == 2)
		{
			if (local_start < nCAP) //bottom
			{
				int size = Min(int(local_n), int(nCAP) - int(local_start));
				for (int i = 0; i < size; i++)
					MPI_Get(&slice[i * n], n, MPI_CXX_DOUBLE_COMPLEX, rank,
							i * n, n, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
			if (local_end > n - 1 - nCAP) //top
			{
				int start = Max(0, (int(n) - int(nCAP)) - int(local_start));
				for (int i = start; i < local_n; i++)
					MPI_Get(&slice[i * n], n, MPI_CXX_DOUBLE_COMPLEX, rank,
							i * n, n, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
		}
		if constexpr (DIM == 3)
		{
			if (local_start < nCAP) //bottom
			{
				int size = Min(int(local_n), int(nCAP) - int(local_start));
				for (int i = 0; i < size; i++)
					for (int j = 0; j < nCAP; j++)
						MPI_Get(&slice[(i * nCAP + j) * n], n, MPI_CXX_DOUBLE_COMPLEX, rank,
								(i * n + (j + boxIndex * nCAP)) * n, n, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
			if (local_end > n - 1 - nCAP) //top
			{
				int start = Max(0, (int(n) - int(nCAP)) - int(local_start));
				for (int i = start; i < local_n; i++)
					for (int j = 0; j < nCAP; j++)
						MPI_Get(&slice[(i * nCAP + j) * n], n, MPI_CXX_DOUBLE_COMPLEX, rank,
								(i * n + (j + boxIndex * nCAP)) * n, n, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
		}
	}
	inline void getYBox(int rank, int boxIndex)
	{
		if constexpr (DIM == 2)
		{
			for (int i = 0; i < local_n; i++)
			{
				MPI_Get(&slice[i * n], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
						i * n, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
				MPI_Get(&slice[(i + 1) * n - nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
						(i + 1) * n - nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}

		}
		if constexpr (DIM == 3)
		{
			for (int i = 0; i < local_n; i++) //left
				for (int j = 0; j < nCAP; j++)
					MPI_Get(&slice[(i * n + j) * nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
							(i * n + j) * n + boxIndex * nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
					// disp, count


			for (int i = 0; i < local_n; i++) //right
				for (int j = n - nCAP; j < n; j++)
					MPI_Get(&slice[(i * n + j) * nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
							(i * n + j) * n + boxIndex * nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
						   //disp, count
		}
	}
	inline void getZBox(int rank, int boxIndex)
	{
		int z_start;
		for (int i = 0; i < local_n; i++)
			for (int j = 0; j < nCAP; j++)
			{
				//front
				z_start = 0;
				MPI_Get(&slice[(i * nCAP + j) * n + z_start], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
						(i * n + (boxIndex * nCAP + j)) * n + z_start, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
				//back
				z_start = n - nCAP;
				MPI_Get(&slice[(i * nCAP + j) * n + z_start], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
						(i * n + (boxIndex * nCAP + j)) * n + z_start, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
				// disp, count
			}
	}

	inline void resetSlice()
	{
		// logInfo("pID %d region %d slice_m %td", pID, region, slice_m);
		for (ind i = 0; i < slice_m; i++)
			slice[i] = zero;

	}
	inline void boxDispatcher(FREE_COORD fc, int rank, int boxIndex)
	{
		// Same nodes in each region talk to each other
		switch (fc)
		{
		case FREE_COORD::X:
			getXBox(rank, boxIndex);
			break;
		case FREE_COORD::Y:
			getYBox(rank, boxIndex);
			break;
		case FREE_COORD::Z:
			getZBox(rank, boxIndex);
			break;
		default:
			break;
		}
	}

	inline double CAP1(ind i) { return absorber.CAP(xmin + i * dx); }
	inline double invCAP1(ind i) { return 1.0 - absorber.CAP(xmin + i * dx); }
	inline double CAP2(ind i, ind j) { return absorber.CAP(xmin + i * dx, xmin + j * dx); }
	inline double invCAP2(ind i, ind j) { return 1 - absorber.CAP(xmin + i * dx, xmin + j * dx); }
	inline void maskRegion(cxd* psi)
	{
		// constexpr auto op = PosOperator();
		// double x;
		ind x, y, z;
		for (ind i = 0; i < local_n; i++)
		{
			if ((sol.freeCoord & FREE_COORD::X) == FREE_COORD::X) x = 0;
			else x = xmin + dx * (i + local_start);
			if constexpr (DIM == 1)
			{
				psi[i] *= absorber.CAP(x);
			}
			else for (ind j = 0; j < n; j++)
			{
				if ((sol.freeCoord & FREE_COORD::Y) == FREE_COORD::Y) y = 0;
				else y = xmin + dx * j;
				if constexpr (DIM == 2)
				{
					psi[i * n + j] *= absorber.CAP(x, y);
				}
				else if constexpr (DIM == 3)
					for (ind k = 0; k < n; k++)
					{
						if ((sol.freeCoord & FREE_COORD::Z) == FREE_COORD::Z) z = 0;
						else z = xmin + dx * k;
						psi[i * nn + j * n + k] *= absorber.CAP(x, y, z);
					}
			}
		}
	}
	constexpr static inline double onesixth = 1.0 / 6.0;
	constexpr static inline double onehalf = 0.5;//0.5;
	static inline bool corrections = true;
	inline void maskSliceX(int boxIndex)
	{
		if constexpr (DIM == 1)
		{
			for (int i = 0; i < local_n; i++)
				slice[i] *= invCAP1(local_start + i);
		}
		if constexpr (DIM == 2)
		{

			for (int i = 0; i < local_n; i++)
				for (int j = 0; j < n; j++)
				{
					slice[(i * n + j)] *= invCAP1(local_start + i);
				}
		}
		if constexpr (DIM == 3)
		{
			for (int i = 0; i < local_n; i++)
				for (int j = 0; j < nCAP; j++)
					for (int k = 0; k < n; k++)
					{
						slice[(i * nCAP + j) * n + k] *= invCAP1(local_start + i);
						if (MPI::group == 1 && corrections)
							slice[(i * nCAP + j) * n + k] *=
							(onesixth * (invCAP2(j + boxIndex * nCAP, k))
							 + onehalf * (CAP1(j + boxIndex * nCAP) + CAP1(k)));

					}
		}
	}
	inline void maskSliceY(int boxIndex)
	{
		if constexpr (DIM == 2)
		{
			for (int i = 0; i < local_n; i++)
				for (int j = 0; j < n; j++)
					slice[(i * n + j)] *= invCAP1(j);
		}
		if constexpr (DIM == 3)
		{
			for (int i = 0; i < local_n; i++)
				for (int j = 0; j < n; j++)
					for (int k = 0; k < nCAP; k++)
					{
						slice[(i * n + j) * nCAP + k] *= invCAP1(j);
						if (MPI::group == 1 && corrections)
							slice[(i * n + j) * nCAP + k] *=
							(onesixth * (invCAP2(local_start + i, boxIndex * nCAP + k))
							 + onehalf * (CAP1(local_start + i) + CAP1(boxIndex * nCAP + k)));
					}
		}
	}
	inline void maskSliceZ(int boxIndex)
	{
		for (int i = 0; i < local_n; i++)
			for (int j = 0; j < nCAP; j++)
				for (int k = 0; k < n; k++)
				{
					slice[(i * nCAP + j) * n + k] *= invCAP1(k);
					if (MPI::group == 1 && corrections)
						slice[(i * nCAP + j) * n + k] *=
						(onesixth * (invCAP2(local_start + i, boxIndex * nCAP + j))
						 + onehalf * (CAP1(local_start + i) + CAP1(boxIndex * nCAP + j)));
				}

	}
	//Max at 0,0,0, Min at (c,c,c), c=CAP_nodes

	inline void coherentAddition(cxd* psi, FREE_COORD fc, int boxIndex)
	{
		if (fc == FREE_COORD::X)
		{
			maskSliceX(boxIndex);
			fftw_execute(transf_slice[0]);
			if constexpr (DIM == 1)
			{
				for (int i = 0; i < local_n; i++)
					psi[i] += slice[i];
			}
			if constexpr (DIM == 2)
			{
				for (int i = 0; i < local_n; i++)
					for (int j = 0; j < n; j++)
						psi[i * n + j] += slice[i * n + j];
			}
			if constexpr (DIM == 3)
			{
				for (int i = 0; i < local_n; i++)
					for (int j = 0; j < nCAP; j++)
						for (int k = 0; k < n; k++)
							psi[(i * n + (j + boxIndex * nCAP)) * n + k] += slice[(i * nCAP + j) * n + k];
			}
		}
		if (fc == FREE_COORD::Y)
		{
			if constexpr (DIM > 1)
			{
				maskSliceY(boxIndex);
				fftw_execute(transf_slice[1]);
				if constexpr (DIM == 2)
				{
					for (int i = 0; i < local_n; i++)
						for (int j = 0; j < n; j++)
							psi[i * n + j] += slice[i * n + j];
				}
				if constexpr (DIM == 3)
				{
					for (int i = 0; i < local_n; i++)
						for (int j = 0; j < n; j++)
							for (int k = 0; k < nCAP; k++)
								psi[(i * n + j) * n + boxIndex * nCAP + k] += slice[(i * n + j) * nCAP + k];
				}
			}
		}

		if (fc == FREE_COORD::Z)
		{
			if constexpr (DIM == 3)
			{
				maskSliceZ(boxIndex);
				fftw_execute(transf_slice[2]);
				for (int i = 0; i < local_n; i++)
					for (int j = 0; j < nCAP; j++)
						for (int k = 0; k < n; k++)
							psi[(i * n + (boxIndex * nCAP + j)) * n + k] += slice[(i * nCAP + j) * n + k];
			}
		}
	}

	inline void transfer()
	{
		Timings::measure::start("TRANSFER");
		for (int gI = 1; gI < MPI::groupCount; gI++)
		{
			// logInfo("Moving to group %d", gI);
			for (int bI = 0; bI < boxesCount; bI++)
			{
				// logInfo("\tDropping box %d", bI);
				for (int sI = 0; sI < gI; sI++)
				{
					// logInfo("\t\t Looking at source number %d", sI);
					if (MPI::group == gI - 1) MPI_Win_fence(MPI_MODE_NOSTORE & MPI_MODE_NOPUT, moreFreeWin);
					if (MPI::group == gI)
					{
						//TODO: check for null ref? 
						auto& src = MPI::sourceRegions[sI];
						resetSlice();
						MPI_Win_fence(MPI_MODE_NOSTORE & MPI_MODE_NOPUT, lessFreeWin);
						boxDispatcher(src.fc, src.rank, bI);
						MPI_Win_fence(MPI_MODE_NOSTORE & MPI_MODE_NOPUT, lessFreeWin);
						coherentAddition(src.fc, bI);
					}
					if (MPI::group == gI - 1) MPI_Win_fence(MPI_MODE_NOSTORE & MPI_MODE_NOPUT, moreFreeWin);
				}
			}
		}
		if (MPI::group != MPI::groupCount - 1) maskRegion();
		Timings::measure::stop("TRANSFER");
	}
};