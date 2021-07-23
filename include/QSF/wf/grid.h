template <class BaseGrid, uind Components, class MPIStrategy = typename BaseGrid::MPIStrategy, class MPIGrids = typename BaseGrid::MPIGrids>
struct Grid;

template <class BaseGrid, uind Components>
struct Grid<BaseGrid, Components, MPI::Slices, MPI::Single> : BaseGrid
{
	// 	static_assert(std::is_base_of_v< CoordinateSystem, CS>,
// 				  "First GridType template argument must be type derived from CoordinateSystem");

// 	// static_assert(is_base_of_v<MPI::Grid, MR>,
// 				//   "Third GridType template argument must be type derived from MPI::Grid");

// 	static_assert((!std::is_same_v<MPI::Rods, MPIStrategy>) || (DIM >= DIMS::D3), "Rod MPI strategy not supported for dimensionality < 3");
	using BaseGrid::pos;
	using BaseGrid::DIM;
	using BaseGrid::absorber;
	using BaseGrid::m;
	using BaseGrid::n;
	using BaseGrid::inv_m;
	using BaseGrid::inv_nn;
	using BaseGrid::inv_n;
	using BaseGrid::nn;
	using BaseGrid::xmin;
	using BaseGrid::dx;
	using BaseGrid::D;
	using MPISol = MPIGrid<typename BaseGrid::MPIGrids, BaseGrid::D, MPI::Slices>;
	static constexpr ind DIMC = DIM + (Components > 1 ? 1 : 0);

	MPISol sol;//MPI Solution
	cxd* psi_total = nullptr; 	// total ψ on m grid (used for saving)
	cxd* psi = nullptr;      	// local ψ(x,t) on m_l grid
	cxd* psi_copy = nullptr;    // backup ψ(x,t) on m_l grid
	cxd* psi_acc = nullptr;     // accumulated ψ(x,t) on m_l grid

	ind m_l; //m local - total number of nodes per process

	ind n0_l; //n0 local - number of n0 nodes attributed to the process
	ind n0_o; //n0 offset - position of first local node in full array
	ind n0_e; //n0 end - position of last local node in full array

	ind n1_l; //n1 local - number of n1 nodes attributed to the process
	ind n1_o; //n1 offset - position of first local node in full array
	ind n1_e; //n1 end - position of last local node in full array

	ind shape[DIMC];/* Shape as given by input, n0*n1*n2*...*/
	ind shape_l[DIMC];/* MPI local shape: local_n0*n1*n2*...*/
	ind shape_t[DIMC]; /*P-space local MPI shape: local_n1*n0*n2*...
	Different from shape_l[] if n0!=n1 and if FFTW_MPI_TRANSPOSED_OUT is used: */

	ind strides[DIMC];
	ind strides_l[DIMC];
	ind strides_t[DIMC];


	ind reverse_shape[DIMC];
	bool mpiFFTW;
	fftw_plan mpi_plans[2];
	fftw_plan extra_plans[2 * DIM];
	int extra_plans_count;

	const int transpose_f[2] = { FFTW_MPI_TRANSPOSED_OUT, FFTW_MPI_TRANSPOSED_IN };
	// const int transpose_f[2] = { 0, 0 };
	const int for_back[2] = { FFTW_FORWARD, FFTW_BACKWARD };



	void gather()
	{
		_logMPI("my %d", MPI::rID);
		if (!MPI::rID && psi_total == nullptr)
		{
			logALLOC("Allocating memory for psi_total");
			psi_total = new cxd[m];
		}
		logMPI("Gathering " psi_symbol "... to address %p", psi_total);
		MPI_Gather(psi, m_l, MPI_CXX_DOUBLE_COMPLEX, psi_total, m_l, MPI_CXX_DOUBLE_COMPLEX, 0, MPI::rComm);
	}


	cxd& operator[](ind index)
	{
		return psi[index];
	}
	const cxd& operator[](ind index) const
	{
		return psi[index];
	}

	template <uind dim = 0, class I, class... Args>
	constexpr inline auto data_offset(I i, Args... args) noexcept
	{
		if constexpr (DIM - 1 == dim) return i * strides_l[dim];
		else return i * strides_l[dim] + data_offset<dim + 1>(args...);
	}
	template <class... I> cxd& operator()(I... i)
	{
		return psi[data_offset(i...)];
	}
	template <class... I> const cxd& operator()(I... i) const
	{
		return psi[data_offset(i...)];
	}

	Grid(Section& settings) : BaseGrid(settings)
	{
		init();
	}

	Grid(BaseGrid g) : BaseGrid(g)
	{
		init();
	}
	void initHelpers()
	{
		logInfo("Using Operator Split Groups (Multiproduct splitting)");
		psi_copy = (cxd*)fftw_malloc(sizeof(cxd) * m_l);
		psi_acc = (cxd*)fftw_malloc(sizeof(cxd) * m_l);
	}
	void init()
	{
		logInfo("Grid Base init %td", n);

		if (!(n % MPI::rSize == 0)) logWarning("Grid length n (%td) should be divisible by the number of MPI region processes (%d) (FFTW reg)", n, MPI::rSize);

		if constexpr (DIM == 1)
		{
			if (!(m % (MPI::rSize * MPI::rSize) == 0))
				logWarning("Grid length n (%td) should be divisible by the number of MPI region processes squared (%d) (FFTW req)", n, MPI::pSize * MPI::pSize);
		}
		else
		{
			if (!((m / n) % (MPI::rSize) == 0))
				logWarning("Row size (m/n=%td) should be divisible by the number of MPI processes (%d) (FFTW req)", m / n, MPI::rSize);
		}

		fftw_mpi_init();
		n0_l = n / MPI::rSize; //Expected split
		n1_l = n0_l; //TODO: Change this to n1/rSize
		// logTestFatal(n >= Im::n && n >= Src::n, "Grid length n (%td) should be larger than IM grid length Im::n (%td) and previous mode grid length Src::n (%td)", n, Im::n, Src::n);

		strides_l[DIM - 1] = 1;
		shape[DIM - 1] = n;
		reverse_shape[0] = n;
		for (int i = DIM - 2; i >= 0; i--)
		{
			shape[i] = n; //Size of each dimension
			reverse_shape[i] = n;
			strides_l[i] = shape[i + 1] * strides_l[i + 1];
		}

		for (int i = 0; i < DIM; i++) logSETUP("Array shape[%d]=%td", i, shape[i]);

		logSETUP("Setting up " psi_symbol "(x,t) arrays and plans...");
		logSETUP("Grid size n (%td) will split into n0_l (%td) by rSize (%d)", n, n0_l, MPI::rSize);

		if constexpr (DIM == 1)
		{
			/* We need to deal with 1D case seperately, see:
			http://fftw.org/doc/Basic-and-advanced-distribution-interfaces.html#Basic-and-advanced-distribution-interfaces */
			ind local_n2, local_start2; //Used in 1D only
			local_n2 = n0_l;
			m_l = fftw_mpi_local_size_1d(n, MPI::rComm, FFTW_FORWARD, MPI::plan_rigor, &n0_l, &n0_o, &local_n2, &local_start2);
			// m_l = fftw_mpi_local_size_1d(n, MPI::rComm, FFTW_BACKWARD, plan_rigor, &n0_l, &n0_o, &local_n2, &local_start2);
		}
		else m_l = fftw_mpi_local_size(DIM, shape, MPI::rComm, &n0_l, &n0_o);
		psi = (cxd*)fftw_malloc(sizeof(cxd) * m_l);
		// We should use non-MPI FFTW routines when transform doesn't involve X (is a free coord)

		// Transposing the first two dimensions in not all-dim cases would lead to problems
		mpiFFTW = sol.bounded[0];
		bool canLeaveTransposed = DIM > 1 && sol.isMain;

		ind nd[1] = { n };

		/* Main region transforms all directions using MPI FFTW, others only use it to transform X */
		if (mpiFFTW)
		{
			//make forward and backward mpi plans
			for (int i = 0; i < 2; i++)
				mpi_plans[i] =
				fftw_mpi_plan_many_dft(sol.isMain ? DIM : 1, //rank
									   sol.isMain ? shape : nd, //dims
									   sol.isMain ? 1 : m / shape[0], //howmany
									   FFTW_MPI_DEFAULT_BLOCK, //block 
									   FFTW_MPI_DEFAULT_BLOCK, //tblock
									   reinterpret_cast<fftw_complex*>(psi),
									   reinterpret_cast<fftw_complex*>(psi),
									   MPI::rComm, for_back[i], MPI::plan_rigor | (canLeaveTransposed ? transpose_f[i] : 0));
		}
		if (mpiFFTW && canLeaveTransposed)
			logWarning("FFTW will leave wf transposed %d %d", transpose_f[0], transpose_f[1]);
		if (sol.isMain)
			logWarning("FFTW will be performed in one step");
		/* Non-mpi fftw (extra) plans for non-main region (region!=0) transforms are coupled together if involve consecutive directions (case for 1-3D)*/
		for (size_t i = 0; i < DIM; i++)
		{
			shape_l[i] = shape[i];
			shape_t[i] = shape[i];
		}
		shape_l[0] = n0_l;
		shape_t[0] = n1_l;
		n1_o = n0_o; ///TODO: change

		n0_e = n0_o + n0_l - 1;
		n1_e = n1_o + n1_l - 1;
		reverse_shape[DIM - 1] = n0_l;

		int fftw_rank = 0;
		int start_dim = -1;
		extra_plans_count = 0;
		if (!sol.isMain) for (int i = 1; i <= DIM; i++)
		{
			if (i == DIM || !sol.bounded[i])
			{
				if (fftw_rank)
				{
					fftw_iodim64* dims = new fftw_iodim64[fftw_rank];
					for (int j = 0; j < fftw_rank; j++)
					{
						dims[j] = { shape_l[start_dim + j], strides_l[start_dim + j],strides_l[start_dim + j] };
						__logMPI("region %d start_dim %d j %d shape_l[start_dim+j]=%td strides_l[start_dim+j]=%td\n", MPI::region, start_dim, j, shape_l[start_dim + j], strides_l[start_dim + j]);
					}

					//define a loop over all howmany_dims (free) dimensions
					int index = 0;
					int howmany_rank = DIM - fftw_rank;
					fftw_iodim64* howmany_dims = new fftw_iodim64[howmany_rank];
					for (int j = 0; j < DIM;j++)
					{
						if (j < start_dim || (j >= start_dim + fftw_rank && j < DIM))
						{
							howmany_dims[index] = { shape_l[j], strides_l[j], strides_l[j] };
							index++;
							__logMPI("region %d start_dim %d index %d j %d shape_l[j]=%td dist[j]=%td\n", MPI::region, start_dim, index, j, shape_l[j], strides_l[j]);
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

		_logMPI("region %d mpiFFTW %d extra_plans_count %d", MPI::region, mpiFFTW, extra_plans_count);




		_logMPI("Region %d MPI::pID: %d got %td nodes or %td rows [start row: %td end row: %td]", MPI::region, MPI::pID, m_l, n0_l, n0_o, n0_e);
		logSETUP(psi_symbol " and forward/backward plans initialized");

		reset();

		// fftw_print_plan(transf_x2p);
		// fftw_print_plan(transf_p2x);
	}
		//Designed to work in X rep, FFTW_MPI_TRANSPOSED_OUT not taken into account
	void add(const cxd* psi_source, ind inp_n, double weight)
	{
		if (n == inp_n)
			for (ind i = 0; i < m_l; i++) psi[i] += weight * psi_source[i];
		else
		{
			ind m_down = (n / 2 - inp_n / 2);
			ind m_up = m_down + inp_n - 1;
			for (ind i = m_down; i <= m_up; i++)
			{
				if (n0_o <= i && i <= n0_e)
				{
					if constexpr (DIM == 1)
						psi[i - n0_o] += weight * psi_source[i - m_down];
					else
					{
						for (ind j = m_down; j <= m_up; j++)
						{
							if constexpr (DIM == 2)
								psi[(i - n0_o) * n + j] += weight * psi_source[(i - m_down) * inp_n + (j - m_down)];
							else
							{
								for (ind k = m_down; k <= m_up; k++)
								{
									psi[(i - n0_o) * nn + j * n + k] += weight * psi_source[((i - m_down) * inp_n + (j - m_down)) * inp_n + (k - m_down)];
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

	inline void backup()
	{
		for (ind i = 0; i < m_l; i++)
		{
			psi_copy[i] = psi[i];
			psi_acc[i] = 0.0;
		}
	}

	inline void restore()
	{
		copyArray(psi, psi_copy);
		// for (ind i = 0; i < m_l; i++) psi[i] = psi_copy[i];
	}

	void accumulate(double coeff)
	{
		for (ind i = 0; i < m_l; i++) psi_acc[i] += coeff * psi[i];
	}

	void collect(double coeff)
	{
		for (ind i = 0; i < m_l; i++)
			psi[i] = coeff * psi[i] + psi_acc[i];
	}

	inline void normalizeAfterTwoFFT()
	{
		switch (sol.freeCoord)
		{
		case FREE_COORD::NO:
			multiplyArray(psi, inv_m); break;
		case FREE_COORD::X:
		case FREE_COORD::Y:
		case FREE_COORD::Z:
			multiplyArray(psi, inv_nn); break;
		case FREE_COORD::XY:
		case FREE_COORD::YZ:
		case FREE_COORD::XZ:
			multiplyArray(psi, inv_n); break;
		default: break;
		}
	}
	int count = 0;
	template <REP R>
	inline void fourier()
	{
		if (count < 2)
		{
			MPI_Barrier(MPI_COMM_WORLD);
			// fprintf(stderr, "REP %s\n\n", R == REP::X ? "P" : "X");
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Barrier(MPI_COMM_WORLD);
				// if (!MPI::rID)
				// {
				// 	gather();
				// 	sprintf(file_path, "%s%d_cpu%d_%s.txt",
				// 			target_path, count, MPI::rSize, R == REP::X ? "P" : "X");
				// 	FILE* f = fopen_with_check<IO_ATTR::WRITE>(file_path);
				// 	for (ind i = 0;i < n;i++)
				// 		for (ind j = 0;j < n;j++)
				// 			fprintf(f, "%td, %td, %g\n", i, j, std::norm(psi_total[i * n + j]));
				// 			// printf("%td, %td, %g\n", i, j, std::norm(psi_total[i * n + j]));
				// 	closeFile(f);
				// }


				// fprintf(stderr, "REP %s\n\n", R == REP::X ? "P" : "X");
			MPI_Barrier(MPI_COMM_WORLD);
			// if (REP::P == R)
			for (int p = 0; p < MPI::rSize;p++)
			{
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
				if (p == MPI::rID)
				{
					if (REP::P == R)
						for (ind i = 0;i < n0_l;i++)
						{
							for (ind j = 0;j < n;j++)
								fprintf(stderr, "%14.2g", 100.0 * std::norm(psi[i * n + j]));
							fprintf(stderr, " <--- %d %s %td\n", MPI::rID, R == REP::X ? "P" : "X", i);
						}

					if (REP::X == R)
						for (ind j = 0;j < n0_l;j++)
						{
							for (ind i = 0;i < n;i++)
								fprintf(stderr, "%14.2d", int(10000.0 * std::norm(psi[j * n + i])));
							fprintf(stderr, " <--- %d %s %td\n", MPI::rID, R == REP::X ? "P" : "X", j);
						}
						// for (ind j = 0;j < n0_l;j++)
						// {
						// 	for (ind i = 0;i < n;i++)
						// 		fprintf(stderr, "%14.2d", int(10000.0 * std::norm(psi[j * n + i])));
						// 	fprintf(stderr, " <--- %d %s %td\n", MPI::rID, R == REP::X ? "P" : "X", j);
						// }
				}
				MPI_Barrier(MPI_COMM_WORLD);
				MPI_Barrier(MPI_COMM_WORLD);
			}

		}
		count++;


	// Timings::measure::start("FFTW");
		static_assert(R == REP::X || R == REP::P,
					  "Can only transform to X or P, unambigously.");
		// logInfo("FFTW into %d", int(R));
		constexpr uind back = R == REP::P ? 0 : 1;
		if (mpiFFTW)
		{
			fftw_execute(mpi_plans[back]);
			// logWarning("mainFFTW");
		}
		for (int i = 0; i < extra_plans_count; i++)
		{
			fftw_execute(extra_plans[i + DIM * back]);
			logWarning("extra_plans");
		}
		if constexpr (back) normalizeAfterTwoFFT();



	// multiplyArray(psi, inv_m * (pow(n, DIM - sol.boundedCoordDim)));
	// Timings::measure::stop("FFTW");
	}

	void post_step()
	{

	}
	void save()
	{
		gather();
		if (!MPI::rID)
		{

			logDUMPS("Dumping " psi_symbol " in X rep");

			FILE* file = openPsi< AFTER<>, REP::X, D, IO_ATTR::WRITE>("name", 0, 0, true);
			//HACK: change second false to true after done with Dmitry!
			writePsiBinaryHeader<D>(file, xmin, -xmin, dx, DUMP_FORMAT{ true, true, true, true, false });
			fwrite(psi_total, sizeof(cxd), m, file);
			fclose(file);
		}
	}
	//BUNCH OF HELPFUL FUNCTIONS
	inline void resetArray(cxd* array)
	{
		for (ind i = 0; i < m_l; i++)
		{
			array[i] = 0.0;
		}
	}
	inline void multiplyArray(cxd* array, const double mult)
	{
		// logInfo("Multiplying array by %g", mult);
		for (ind i = 0; i < m_l; i++)
			array[i] *= mult;
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
};




template <class BaseGrid, uind Components>
struct Grid<BaseGrid, Components, MPI::Slices, MPI::Multi> : Grid<BaseGrid, Components, MPI::Slices, MPI::Single>
{
	using Base = Grid<BaseGrid, Components, MPI::Slices, MPI::Single>;
	using Base::DIM;
	using Base::pos;
	using Base::sol;
	using Base::psi;
	using Base::m;
	using Base::n;
	using Base::xmin;
	using Base::dx;
	using Base::nn;
	using Base::n0_l;
	using Base::absorber;
	using Base::n0_o;
	using Base::n0_e;
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
	void init()
	{
		logInfo("Grid MultiExtension init");
		// absorber.correct(n, L);
		// nCAP = absorber.nCAP;
		// boxesCount = DIM < 3 ? 1 : n / nCAP;
		// CAPnodesPerP = n0_l / nCAP;
		initSlice();
		initWindow();
	}
	Grid(Section& settings) :Base(settings) { init(); }

	Grid(BaseGrid g) : Base(g) { init(); }

	// ~Grid()
	// {
	// 	// destroySlice();
	// }

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

		ind shape[DIM];
		shape[0] = n0_l;
		ind strides_l[DIM];
		strides_l[DIM - 1] = 1;
		for (int i = 1; i < DIM; i++)
		{
			//HACK: very dirty
			if constexpr (DIM == 2) shape[1] = n;
			else if constexpr (DIM == 3)
			{
				shape[1] = (i == 1 ? n : nCAP);
				shape[2] = (i == 1 ? nCAP : n);
			}

			for (int j = DIM - 2; j >= 0; j--)
				strides_l[j] = shape[j + 1] * strides_l[j + 1];

			fftw_iodim64 dims[1] = { n , strides_l[i], strides_l[i] };
			fftw_iodim64 howmany_dims[DIM - 1];
			int index = 0;
			for (int j = 0; j < DIM; j++)
			{
				if (j != i)
				{
					howmany_dims[index] = { shape[j], strides_l[j],strides_l[j] };
					index++;
				}
			}
			printf("i=%d shape[i]=%td %td\n", i, shape[i], strides_l[i]);
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
			if (n0_o < nCAP) //bottom
			{
				int size = Min(int(n0_l), int(nCAP) - int(n0_o));
				MPI_Get(&slice[0], size, MPI_CXX_DOUBLE_COMPLEX, rank,
						0, size, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
			if (n0_e > n - 1 - nCAP) //top
			{
				int start = Max(0, int(n) - int(nCAP) - int(n0_o));
				MPI_Get(&slice[start], n0_l - start, MPI_CXX_DOUBLE_COMPLEX, rank,
						start, n0_l - start, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
		}
		if constexpr (DIM == 2)
		{
			if (n0_o < nCAP) //bottom
			{
				int size = Min(int(n0_l), int(nCAP) - int(n0_o));
				for (int i = 0; i < size; i++)
					MPI_Get(&slice[i * n], n, MPI_CXX_DOUBLE_COMPLEX, rank,
							i * n, n, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
			if (n0_e > n - 1 - nCAP) //top
			{
				int start = Max(0, (int(n) - int(nCAP)) - int(n0_o));
				for (int i = start; i < n0_l; i++)
					MPI_Get(&slice[i * n], n, MPI_CXX_DOUBLE_COMPLEX, rank,
							i * n, n, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
		}
		if constexpr (DIM == 3)
		{
			if (n0_o < nCAP) //bottom
			{
				int size = Min(int(n0_l), int(nCAP) - int(n0_o));
				for (int i = 0; i < size; i++)
					for (int j = 0; j < nCAP; j++)
						MPI_Get(&slice[(i * nCAP + j) * n], n, MPI_CXX_DOUBLE_COMPLEX, rank,
								(i * n + (j + boxIndex * nCAP)) * n, n, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}
			if (n0_e > n - 1 - nCAP) //top
			{
				int start = Max(0, (int(n) - int(nCAP)) - int(n0_o));
				for (int i = start; i < n0_l; i++)
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
			for (int i = 0; i < n0_l; i++)
			{
				MPI_Get(&slice[i * n], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
						i * n, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
				MPI_Get(&slice[(i + 1) * n - nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
						(i + 1) * n - nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
			}

		}
		if constexpr (DIM == 3)
		{
			for (int i = 0; i < n0_l; i++) //left
				for (int j = 0; j < nCAP; j++)
					MPI_Get(&slice[(i * n + j) * nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
							(i * n + j) * n + boxIndex * nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
					// disp, count


			for (int i = 0; i < n0_l; i++) //right
				for (int j = n - nCAP; j < n; j++)
					MPI_Get(&slice[(i * n + j) * nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
							(i * n + j) * n + boxIndex * nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
						   //disp, count
		}
	}
	inline void getZBox(int rank, int boxIndex)
	{
		int z_start;
		for (int i = 0; i < n0_l; i++)
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

	inline double CAP1(ind i) { return absorber(xmin + i * dx); }
	inline double invCAP1(ind i) { return 1.0 - absorber(xmin + i * dx); }
	inline double CAP2(ind i, ind j) { return absorber.CAP(xmin + i * dx, xmin + j * dx); }
	inline double invCAP2(ind i, ind j) { return 1 - absorber(xmin + i * dx, xmin + j * dx); }
	inline void maskRegion()
	{
		// constexpr auto op = PosOperator();
		// double x;
		ind x, y, z;
		for (ind i = 0; i < n0_l; i++)
		{
			if ((sol.freeCoord & FREE_COORD::X) == FREE_COORD::X) x = 0;
			else x = xmin + dx * (i + n0_o);
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
						psi[i * nn + j * n + k] *= absorber(x, y, z);
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
			for (int i = 0; i < n0_l; i++)
				slice[i] *= invCAP1(n0_o + i);
		}
		if constexpr (DIM == 2)
		{

			for (int i = 0; i < n0_l; i++)
				for (int j = 0; j < n; j++)
				{
					slice[(i * n + j)] *= invCAP1(n0_o + i);
				}
		}
		if constexpr (DIM == 3)
		{
			for (int i = 0; i < n0_l; i++)
				for (int j = 0; j < nCAP; j++)
					for (int k = 0; k < n; k++)
					{
						slice[(i * nCAP + j) * n + k] *= invCAP1(n0_o + i);
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
			for (int i = 0; i < n0_l; i++)
				for (int j = 0; j < n; j++)
					slice[(i * n + j)] *= invCAP1(j);
		}
		if constexpr (DIM == 3)
		{
			for (int i = 0; i < n0_l; i++)
				for (int j = 0; j < n; j++)
					for (int k = 0; k < nCAP; k++)
					{
						slice[(i * n + j) * nCAP + k] *= invCAP1(j);
						if (MPI::group == 1 && corrections)
							slice[(i * n + j) * nCAP + k] *=
							(onesixth * (invCAP2(n0_o + i, boxIndex * nCAP + k))
							 + onehalf * (CAP1(n0_o + i) + CAP1(boxIndex * nCAP + k)));
					}
		}
	}
	inline void maskSliceZ(int boxIndex)
	{
		for (int i = 0; i < n0_l; i++)
			for (int j = 0; j < nCAP; j++)
				for (int k = 0; k < n; k++)
				{
					slice[(i * nCAP + j) * n + k] *= invCAP1(k);
					if (MPI::group == 1 && corrections)
						slice[(i * nCAP + j) * n + k] *=
						(onesixth * (invCAP2(n0_o + i, boxIndex * nCAP + j))
						 + onehalf * (CAP1(n0_o + i) + CAP1(boxIndex * nCAP + j)));
				}

	}
	//Max at 0,0,0, Min at (c,c,c), c=CAP_nodes

	inline void coherentAddition(FREE_COORD fc, int boxIndex)
	{
		if (fc == FREE_COORD::X)
		{
			maskSliceX(boxIndex);
			fftw_execute(transf_slice[0]);
			if constexpr (DIM == 1)
			{
				for (int i = 0; i < n0_l; i++)
					psi[i] += slice[i];
			}
			if constexpr (DIM == 2)
			{
				for (int i = 0; i < n0_l; i++)
					for (int j = 0; j < n; j++)
						psi[i * n + j] += slice[i * n + j];
			}
			if constexpr (DIM == 3)
			{
				for (int i = 0; i < n0_l; i++)
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
					for (int i = 0; i < n0_l; i++)
						for (int j = 0; j < n; j++)
							psi[i * n + j] += slice[i * n + j];
				}
				if constexpr (DIM == 3)
				{
					for (int i = 0; i < n0_l; i++)
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
				for (int i = 0; i < n0_l; i++)
					for (int j = 0; j < nCAP; j++)
						for (int k = 0; k < n; k++)
							psi[(i * n + (boxIndex * nCAP + j)) * n + k] += slice[(i * nCAP + j) * n + k];
			}
		}
	}

	void transfer()
	{
		// Timings::measure::start("TRANSFER");
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
		// Timings::measure::stop("TRANSFER");
	}

	void post_step()
	{
		Base::post_step();
		transfer();
	}
};