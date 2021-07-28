
template <class Hamiltonian, class BaseGrid, uind Components, class MPI_GC>
struct LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, true> :
	LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, false>
{
	using Base = LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, false>;
	using Base::DIM, Base::n_lx;
	// using Base::n_lx;
	using Base::m_l;
	using Base::mcomm;
	using Base::psi;
	using Base::BaseGrid;

	static inline constexpr cxd zero = { 0.0, 0.0 };
	ind slice_n;
	ind slice_start;
	ind slice_n2;
	ind slice_start2;
	ind slice_m;
	fftw_plan transf_slice[DIM];

	MPI_Win lessFreeWin = nullptr;
	MPI_Win moreFreeWin = nullptr;

	cxd* slice;
	ind sliceSize;
	int boxesCount;
	int nCAP;
	bool CAPnodesPerP;

	static_assert(Base::hasAbsorber, "Absorber is required for multigrid computations");
	LocalGrid(Base g) : Base(g) { init(); }
	LocalGrid(Section& settings) :Base(settings) { init(); }

	void init()
	{
		logInfo("MultiRegions init");
		// absorber.correct(n, L);
		// nCAP = absorber.nCAP;
		// boxesCount = DIM < 3 ? 1 : n / nCAP;
		// CAPnodesPerP = n0_l / nCAP;
		initSlice();
		if (mcomm.moreFree) MPI_Win_create(psi, m_l * sizeof(cxd), sizeof(cxd), MPI_INFO_NULL, mcomm.moreFree, &moreFreeWin);
		if (mcomm.lessFree) MPI_Win_create(NULL, 0, sizeof(cxd), MPI_INFO_NULL, mcomm.lessFree, &lessFreeWin);
	}

	~LocalGrid()
	{
		if (!mcomm.isMain) delete[] slice;
		if (mcomm.moreFree) MPI_Win_free(&moreFreeWin);
		if (mcomm.lessFree) MPI_Win_free(&lessFreeWin);
		~Base::LocalGrid();
	}

	void initSlice()
	{
		/* TODO: Here we decide how the wf will be sliced.
		1) SIZE: In 1D slices are not cut, hence sliceSize == n.
		In other cases size of the slice should be constant for all regions sliceSize == nCAP*m/n.
		2) Here we want to split the job among all the processes in the region, hence the slices
		will be cut along the X-direction, whenever we are transfering a direction != X.
		When transfering X direction we cut along Y direction as we need to perform fftw along X. */
		if (mcomm.isMain) return; //Main region doesn't operate with slices

		sliceSize = DIM == 1 ? n_lx[0] : (DIM == 2 ? n_lx[0] * n_lx[1] : nCAP * n_lx[0] * n_lx[1]);
		ind howmany = DIM == 1 ? 1 : (DIM == 2 ? n : nCAP * n);
		ind dims[1] = { n };
		slice = (cxd*)fftw_malloc(sizeof(cxd) * sliceSize / MPI::rSize);
		transf_slice[0] = fftw_mpi_plan_many_dft(1, dims, howmany,
												 FFTW_MPI_DEFAULT_BLOCK, //block 
												 FFTW_MPI_DEFAULT_BLOCK, //tblock
												 reinterpret_cast<fftw_complex*>(slice),
												 reinterpret_cast<fftw_complex*>(slice),
												 MPI::rComm, FFTW_FORWARD, MPI::plan_rigor);
		if (transf_slice[0] == NULL) logInfo("fftw_mpi_plan... (slice) returned NULL");

		ind n_l[DIM];
		n_l[0] = n0_l;
		ind strides_l[DIM];
		strides_l[DIM - 1] = 1;
		for (int i = 1; i < DIM; i++)
		{
			//HACK: very dirty
			if constexpr (DIM == 2) n_l[1] = n;
			else if constexpr (DIM == 3)
			{
				n_l[1] = (i == 1 ? n : nCAP);
				n_l[2] = (i == 1 ? nCAP : n);
			}

			for (int j = DIM - 2; j >= 0; j--)
				strides_l[j] = n_l[j + 1] * strides_l[j + 1];

			fftw_iodim64 dims[1] = { n , strides_l[i], strides_l[i] };
			fftw_iodim64 howmany_dims[DIM - 1];
			int index = 0;
			for (int j = 0; j < DIM; j++)
			{
				if (j != i)
				{
					howmany_dims[index] = { n_l[j], strides_l[j],strides_l[j] };
					index++;
				}
			}
			printf("i=%d n_l[i]=%td %td\n", i, n_l[i], strides_l[i]);
			transf_slice[i] = fftw_plan_guru64_dft(
				1, dims, DIM - 1, howmany_dims, //rank, dims, howmany_rank, howmany_dims
				reinterpret_cast<fftw_complex*>(slice),
				reinterpret_cast<fftw_complex*>(slice),
				FFTW_FORWARD, MPI::plan_rigor);
		}
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

#pragma region Tests
		//See: http://www.fftw.org/fftw3_doc/Load-balancing.html#Load-balancing
	void test()
	{
		// logTestFatal(n % MPI::rSize == 0, "Grid length n (%td) should be divisible by the number of MPI region processes (%d) (FFTW reg)", n, MPI::rSize);
		if constexpr (DIM == 1)
		{
			// logTestFatal(m % (MPI::rSize * MPI::rSize) == 0, "Grid length n (%td) should be divisible by the number of MPI region processes squared (%d) (FFTW req)", n, MPI::pSize * MPI::pSize);
		}
		else
		{
			// logTestFatal((m / n) % (MPI::rSize) == 0, "Row size (m/n=%td) should be divisible by the number of MPI processes (%d) (FFTW req)", m / n, MPI::rSize);
		}

		// logTestFatal(n >= Im::n && n >= Src::n, "Grid length n (%td) should be larger than IM grid length Im::n (%td) and previous mode grid length Src::n (%td)", n, Im::n, Src::n);
	}
#pragma endregion Tests

};