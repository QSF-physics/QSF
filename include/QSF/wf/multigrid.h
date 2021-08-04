struct SrcRegion
{
	int rank;
	AXIS fc;
};
template <class Hamiltonian, class BaseGrid, uind Components, class MPI_GC>
struct LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, true> :
	LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, false>
{
	using Base = LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, false>;
	using Base::DIM, Base::n_lx, Base::n, Base::pos_lx, Base::nCAP;
	// using Base::n_lx;
	using Base::m_l, Base::m, Base::xmin, Base::dx;
	using Base::mcomm, Base::psi, Base::mask;
	static_assert(Base::hasAbsorber, "Absorb is required for multigrid computations");
	static inline constexpr cxd zero = { 0.0, 0.0 };

	fftw_plan transf_slice[DIM];
	MPI_Win lessFreeWin = nullptr;
	MPI_Win moreFreeWin = nullptr;

	cxd* slice;
	ind m_slice;
	ind m_lslice;
	int boxesCount;

	bool isBottomX;
	bool isTopX;
	ind sizeX;
	ind startX;
	ind stride_slice[DIM][DIM]; //For each "free" axis
	template <ind Is> ind static constexpr rev = DIM - 2 - Is;

	std::vector<SrcRegion> sourceRegions;

	LocalGrid(Base g) : Base(g) { init(); }
	LocalGrid(Section& settings) :Base(settings) { init(); }

	void init()
	{
		logInfo("MultiRegions init");
		int index = 0;
		for (int i = 0; i < MPI::regionCount; i++)
		{
			if (mcomm.freeCoordsCount[MPI::region] - mcomm.freeCoordsCount[i] == 1)
			{
				switch (mcomm.freeCoord ^ mcomm.freeCoords[i])
				{
				case AXIS::X: sourceRegions.push_back({ index * MPI::rSize + MPI::rID,AXIS::X }); break;
				case AXIS::Y: sourceRegions.push_back({ index * MPI::rSize + MPI::rID,AXIS::Y }); break;
				case AXIS::Z: sourceRegions.push_back({ index * MPI::rSize + MPI::rID,AXIS::Z }); break;
				default: break;
				}
				index++;
			}
		}

		_logMPI("[group %d region %d pID %d] has %td sources", MPI::group, MPI::region, MPI::pID, sourceRegions.size());

/* 		if (n % nCAP)
		{
			ind lCAP = nCAP - 1;
			ind rCAP = nCAP + 1;
			while ((n % lCAP) && (n % rCAP)) { lCAP--; rCAP++; }
			if ((rCAP > n2[0] - 1) && (lCAP < 1))
			{
				logError("Auto enlarging and decreasing nCAP failed.");
			}
			else if (rCAP > n2 - 1)
			{
				logWarning("Auto enlarging nCAP failed.");
			}
			else if (lCAP < 1)
			{
				logWarning("Auto decreasing nCAP failed.");
			}
			else
			{
				if (n % rCAP == 0)
				{
					nCAP = rCAP;
					CAPlength = nCAP * dx;
					logWarning("Auto enlarging nCAP to %td, new CAPlength %g", nCAP, CAPlength);
				}
				else if (n % lCAP == 0)
				{
					nCAP = lCAP;
					CAPlength = nCAP * dx;
					logWarning("Auto decreasing nCAP to %td, new CAPlength %g", nCAP, CAPlength);
				}
			}
		} */
		boxesCount = DIM < 3 ? 1 : n[0] / nCAP;

		logWarning("boxesCount %d", boxesCount);
		if (!mcomm.isMain) initSlice(n_seq<DIM>); //Main region doesn't operate with slices
		if (mcomm.moreFree) MPI_Win_create(psi, m_l * sizeof(cxd), sizeof(cxd), MPI_INFO_NULL, mcomm.moreFree, &moreFreeWin);
		if (mcomm.lessFree) MPI_Win_create(NULL, 0, sizeof(cxd), MPI_INFO_NULL, mcomm.lessFree, &lessFreeWin);
	}

	~LocalGrid()
	{
		logWarning("localGrid destructor called");
		// if (!mcomm.isMain) delete[] slice;
		// MPI_Barrier(MPI_COMM_WORLD);
		// if (mcomm.moreFree) MPI_Win_free(&moreFreeWin);
		// if (mcomm.lessFree) MPI_Win_free(&lessFreeWin);
		// ~Base();
	}

	template <uind dirFree, uind ...dirs>
	void initSliceStrides(seq<dirs...>)
	{
		stride_slice[dirFree][DIM - 1] = 1;
		((stride_slice[dirFree][dirs] = sliceShape<dirFree, dirs + 1>() * stride_slice[dirFree][dirs + 1]), ...);
		std::cout << MPI::pID << " dirFree: " << dirFree << " stride_slice[dirFree][..]" << stride_slice[dirFree][0] << " " << stride_slice[dirFree][1] << " " << stride_slice[dirFree][2] << " " << std::endl;
	}
	// Here we initialize the slices (temp wf containers used for wf transfering)
	// By design we want them to be universal, i.e. hold pieces of the wf coming 
	// from sources with one less free coordinate and apply FFTW to them, hence
	// size is (before MPI splitting) in 1D: m=n[0], 2D: m=n[0]*n[1], 3D: nCAP*m/n
	// or 1D: m_l=n_lx[0], 2D: m_l=n_lx[0]*n[1], 3D: nCAP*n_lx[0]*n[0], where 
	// n_lx[0]=n[0]/MPI::rSize. This means that in the 1D and 2D case we in fact 
	// transfer all the wf at once. It would be possible to save some space in 2D, 
	// however the per-process access to the wf pieces would more complicated.
	// This implementation is not the most efficient, but it splits the job among 
	// all the processes in the region. In 3D case the nCAP cuts are NEVER made 
	// along the X dir (which is split by fftw_mpi anyway).
	template <uind dirFree, uind ...dirs>
	void initFFTWplans(seq<dirs...>)
	{

		if constexpr (dirFree == 0)
		{
			printf("pID %d MPI initFFTWplans dirFree: [%td] howmany %td\n",
				   MPI::pID, dirFree, m_slice / n[0]);
			transf_slice[0] =
				fftw_mpi_plan_many_dft(1, n, m_slice / n[0],
									   FFTW_MPI_DEFAULT_BLOCK, //block 
									   FFTW_MPI_DEFAULT_BLOCK, //tblock
									   reinterpret_cast<fftw_complex*>(slice),
									   reinterpret_cast<fftw_complex*>(slice),
									   MPI::rComm, FFTW_FORWARD, MPI::plan_rigor);
		}
		else
		{
			printf("pID %d initFFTWplans dirFree: [%td] howmany %td\n", MPI::pID, dirFree, m_slice / n[0]);
			// The dimension we're going to transform
			fftw_iodim64 dims[1] = { {n[dirFree] , stride_slice[dirFree][dirFree], stride_slice[dirFree][dirFree]} };
			fftw_iodim64 howmany_dims[DIM - 1];
			int index = 0;
			if (MPI::pID == 1)
				std::cout << "---:> pID:" << MPI::pID << " dirFree:" << dirFree << " dims: " << dims[0].n << " " << dims[0].is << std::endl;

			([&] {if (dirs != dirFree) {
				howmany_dims[index] = { sliceShape<dirFree, dirs>(), stride_slice[dirFree][dirs], stride_slice[dirFree][dirs] };
				index++;
				if (MPI::pID == 1)
					std::cout << "---> pID:" << MPI::pID << " dirFree:" << dirFree << " dirs:" << dirs << "howmany_dims: " << howmany_dims[index - 1].n << " " << howmany_dims[index - 1].is << std::endl;
			}}(), ...);

			transf_slice[dirFree] = fftw_plan_guru64_dft(
				1, dims, DIM - 1, howmany_dims, //rank, dims, howmany_rank, howmany_dims
				reinterpret_cast<fftw_complex*>(slice),
				reinterpret_cast<fftw_complex*>(slice),
				FFTW_FORWARD, MPI::plan_rigor);
		}
	}


	template<uind ... dirFree>
	void initSlice(seq<dirFree...>)
	{
		m_slice = DIM == 3 ? m / n[2] * nCAP : m;
		m_lslice = m_slice / MPI::rSize;
		slice = (cxd*)fftw_malloc(sizeof(cxd) * m_lslice);
		printf("%d initSlice m_slice [%td] m_lslice [%td]\n", MPI::pID, m_slice, m_lslice);

		(initSliceStrides<dirFree>(rev_seq<DIM - 1>), ...);
		(initFFTWplans<dirFree>(n_seq<DIM>), ...);

		isBottomX = pos_lx.first < nCAP;
		sizeX = Min(n_lx[0], nCAP - pos_lx.first);

		isTopX = pos_lx.last > n[0] - 1 - nCAP;
		startX = Max(0, n[0] - nCAP - pos_lx.first);

		printf("pID: %d  startX[%td] sizeX[%td]\n", MPI::pID, startX, sizeX);
		// if (transf_slice[0] == NULL) logInfo("fftw_mpi_plan... (slice) returned NULL");
	}

	template <uind dirFree, uind dir>
	inline bool slidesIn()
	{
		if constexpr (DIM == 3)
		{
			if constexpr (dir == 1 && (Axis<dirFree> == AXIS::X || Axis<dirFree> == AXIS::Z))
				return true;
			if constexpr (dir == 2 && Axis<dirFree> == AXIS::Y)
				return true;
		}
		return false;
	}

	template <uind dirFree, uind dir>
	inline constexpr ind sliceShape()
	{
		//transfer X/Z ==> box shift over Y, transfer Y ==> box shift over Z
		if (slidesIn<dirFree, dir>()) return nCAP;
		else return Base::template shape<REP::X, dir>();
	}

	template <uind dirFree, uind dim = 0, class I, class... Args>
	constexpr inline auto slice_offset(I i, Args... args) const// noexcept
	{
		if constexpr (DIM - 1 == dim) return i * stride_slice[dirFree][dim];
		else return i * stride_slice[dirFree][dim] + slice_offset<dirFree, dim + 1>(args...);
	}

	template <uind dirFree, uind dir>
	inline ind slider(ind boxIndex)
	{
		return slidesIn<dirFree, dir>() ? boxIndex * nCAP : 0;
	}

	template <uind dirFree, uind dir>
	inline constexpr ind sliceLimited()
	{
		if constexpr (dirFree == dir && dirFree) return nCAP;
		else return sliceShape<dirFree, dir>();
	}
	template <uind dirFree, uind dir>
	inline ind top()
	{
		return ((dirFree == dir && dirFree) ? n_lx[dir] - nCAP : 0);
	}

	template <uind dirFree, uind dir>
	inline ind tslider(ind boxIndex)
	{
		return slider<dirFree, dir>(boxIndex) + top<dirFree, dir>();
	}
	inline void getData(int rank, ind slice_off, ind wf_offset, ind count)
	{
		MPI_Get(&slice[slice_off], //pos in slice
				count, //send count
				MPI_CXX_DOUBLE_COMPLEX, rank,
				wf_offset, //pos in wf
				count, //recv count
				MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
	}


	template <uind dirFree, uind ... dirs>
	inline void getBox(int rank, int box, seq<dirs...>)
	{
		// printf("TTT %td %td\n", sliceLimited<dirFree, rev<dirs>>()...);
		ind counters[DIM] = { 0 }; //Z dimension goes in "bulk" MPI_Get
		do
		{

			// negative axis
			if (Axis<dirFree> != AXIS::X || (isBottomX && counters[0] < sizeX))
				getData(rank, slice_offset<dirFree>(counters[dirs]..., 0),
						Base::template data_offset<REP::X>(
							(counters[dirs] + slider<dirFree, dirs>(box))..., slider<dirFree, DIM - 1>(box)),
						sliceLimited<dirFree, DIM - 1>());

			// positive axis
			if (Axis<dirFree> != AXIS::X || (isTopX && counters[0] >= startX))
				getData(rank,
						slice_offset<dirFree>((counters[dirs] + top<dirFree, dirs>())..., top<dirFree, DIM - 1>()),
						Base::template data_offset<REP::X>((counters[dirs] + tslider<dirFree, dirs>(box))..., tslider<dirFree, DIM - 1>(box)),
						sliceLimited<dirFree, DIM - 1>());
			// if ((Axis<dirFree> != AXIS::X || (isTopX && counters[0] >= startX)) && MPI::pID == 1)
				// printf("==== %td %td %td\n", counters[dirs] ..., startX);

		} while (!(...&&
				   ((counters[rev<dirs>]++, counters[rev<dirs>] < sliceLimited<dirFree, rev<dirs>>())
					? (false) : (counters[rev<dirs>] = 0, true))
				   ));
	}

	inline void resetSlice()
	{
		// logInfo("pID %d region %d m_lslice %td", pID, region, m_lslice);
		for (ind i = 0; i < m_lslice; i++) slice[i] = zero;
	}
	inline void boxDispatcher(AXIS fc, int rank, int boxIndex)
	{
		// std::cout << " -->box dispatch" << rank << " " << boxIndex << " fc:" << uind(fc) << " pID" << MPI::pID << std::endl;
		// Same nodes in each region talk to each other
		// getXBox(renk)
		switch (fc)
		{
		case AXIS::X:
			getBox<0>(rank, boxIndex, n_seq<DIM - 1>);
			// getXBox(rank, boxIndex);
			break;
		case AXIS::Y:
			getBox<1>(rank, boxIndex, n_seq<DIM - 1>);
			// getYBox(rank, boxIndex);
			break;
		case AXIS::Z:
			getBox<2>(rank, boxIndex, n_seq<DIM - 1>);
			// getZBox(rank, boxIndex);
			break;
		default:
			break;
		}
	}


	inline double CAP1(ind i) { return Base::template mask<0>(i); }
	inline double invCAP1(ind i) { return 1.0 - Base::template mask<0>(i); }
	inline double CAP2(ind i, ind j) { return Base::template mask<0, 1>(i, j); }
	inline double invCAP2(ind i, ind j) { return 1 - Base::template mask<0, 1>(i, j); }



	inline static constexpr bool corrections = true;

	// template <uind dirFree, uind ... dirs>
	// double weightCorrections(ind x, ind y, ind z)
	// {
	// 	return (onesixth * (invCAP2(pos_lx.first + i, boxIndex * nCAP + k))
	// 			+ onehalf * (CAP1(pos_lx.first + i) + CAP1(boxIndex * nCAP + k)));
	// }
	const double onesixth = 1.0 / 6.0;
	const double onehalf = 0.5;//0.5;
	template <uind dirFree, uind ... dirs, typename ...Nodes>
	inline double mask_correction(Nodes ... nodes)
	{
		return 1.0;
		if (MPI::group == 1) return 1.0;
		else return  onehalf * ((dirs == dirFree ? 0 : Base::template mask<dirs>(nodes)) + ... + 0) +
			onesixth * (1 - ((dirs == dirFree ? 1.0 : Base::template mask<dirs>(nodes))*...));
	}

	template <uind dirFree, uind ... dirs>
	inline void maskSlice(int boxIndex, seq<dirs...>)
	{
		ind counters[DIM] = { 0 };
		do
		{
			//those by design operate on nCAP nodes, which takes neg dist from boundry
			if (Axis<dirFree> != AXIS::X || (isBottomX && counters[0] < sizeX))
				slice[slice_offset<dirFree>(counters[dirs]...)]
				*= Base::template inv_mask<0>(Base::template neg_dist_from_edge<REP::X, dirFree>(counters[dirFree]))
				* mask_correction<dirFree, dirs...>(Base::template neg_dist_from_edge<REP::X, dirFree>(counters[dirs] + slider<dirFree, dirs>(boxIndex))...);

			if (Axis<dirFree> != AXIS::X || (isTopX && counters[0] >= startX))
				slice[slice_offset<dirFree>((counters[dirs] + top<dirFree, dirs>())...)]
				*= Base::template inv_mask<0>(Base::template neg_dist_from_edge<REP::X, dirFree>(counters[dirFree] + top<dirFree, dirFree>()))
				* mask_correction<dirFree, dirs...>(Base::template neg_dist_from_edge<REP::X, dirFree>(counters[dirs] + tslider<dirFree, dirs>(boxIndex))...);


		} while (!(...&&
				   ((counters[Base::template rev<dirs>]++,
					 counters[Base::template rev<dirs>] < sliceLimited<dirFree, Base::template rev<dirs>>())
					? (false) : (counters[Base::template rev<dirs>] = 0, true))
				   ));
	}

	template <uind dirFree, uind ... dirs>
	inline void add(int box, seq<dirs...>)
	{
		ind counters[DIM] = { 0 };
		do
		{
			psi[Base::template data_offset<REP::X>((counters[dirs] + slider<dirFree, dirs>(box))...)]
				+= slice[slice_offset<dirFree>(counters[dirs]...)];

		} while (!(...&&
				   ((counters[Base::template rev<dirs>]++,
					 counters[Base::template rev<dirs>] < sliceShape<dirFree, Base::template rev<dirs>>())
					? (false) : (counters[Base::template rev<dirs>] = 0, true))
				   ));
	}

	//Max at 0,0,0, Min at (c,c,c), c=CAP_nodes
	inline void coherentAddition(AXIS fc, int boxIndex)
	{
		// printf("cohh\n");
		if (fc == AXIS::X)
		{
			maskSlice<0>(boxIndex, n_seq<DIM>);
			fftw_execute(transf_slice[0]);
			add<0>(boxIndex, n_seq<DIM>);
		}
		if (fc == AXIS::Y)
		{
			if constexpr (DIM > 1)
			{
				maskSlice<1>(boxIndex, n_seq<DIM>);
				fftw_execute(transf_slice[1]);
				add<1>(boxIndex, n_seq<DIM>);
			}
		}

		if (fc == AXIS::Z)
		{
			if constexpr (DIM > 2)
			{
				maskSlice<2>(boxIndex, n_seq<DIM>);
				fftw_execute(transf_slice[2]);
				add<2>(boxIndex, n_seq<DIM>);
			}
		}
	}
	void postCompute()
	{
		Base::postCompute();
		transfer();
		if (MPI::group != MPI::groupCount - 1) Base::maskRegion();
	}
	void transfer()
	{
		// MPI_Barrier(MPI_COMM_WORLD);
		// Timings::measure::start("TRANSFER");
		for (int gI = 1; gI < MPI::groupCount; gI++) //groups
		{
			// logInfo("Moving to group %d", gI);
			for (int bI = 0; bI < boxesCount; bI++) //boxes
			{
				// logInfo("\tDropping box %d", bI);
				for (int sI = 0; sI < gI; sI++) //sources
				{
					// logInfo("\t\t Looking at source number %d", sI);
					// std::cout << "pID: " << MPI::pID << " group: " << MPI::group << " region: " << MPI::region << " gI:" << gI << " bI: " << bI << " sI: " << sI << std::endl;
					if (MPI::group == gI - 1) MPI_Win_fence(MPI_MODE_NOSTORE & MPI_MODE_NOPUT, moreFreeWin);
					if (MPI::group == gI)
					{
						//TODO: check for null ref?
						auto& src = sourceRegions[sI];
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




// // (0,0,0) = (bottom, left, front)
// 	inline void getXBox(int rank, int boxIndex)
// 	{
// 		if constexpr (DIM == 1)
// 		{
// 			if (isBottomX) //bottom
// 				MPI_Get(&slice[0], sizeX, MPI_CXX_DOUBLE_COMPLEX, rank,
// 						0, sizeX, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);

// 			if (isTopX) //top
// 				MPI_Get(&slice[startX], n_lx[0] - startX, MPI_CXX_DOUBLE_COMPLEX, rank,
// 						startX, n_lx[0] - startX, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 		}
// 		if constexpr (DIM == 2)
// 		{
// 			if (isBottomX) //bottom
// 				for (int i = 0; i < sizeX; i++)
// 					MPI_Get(&slice[i * n[0]], n[0], MPI_CXX_DOUBLE_COMPLEX, rank,
// 							i * n[0], n[0], MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);

// 			if (isTopX) //top
// 				for (int i = startX; i < n_lx[0]; i++)
// 					MPI_Get(&slice[i * n[0]], n[0], MPI_CXX_DOUBLE_COMPLEX, rank,
// 							i * n[0], n[0], MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 		}
// 		if constexpr (DIM == 3)
// 		{
// 			// std::cout << ("get box") << std::endl;
// 			if (isBottomX) //bottom
// 				for (int i = 0; i < sizeX; i++) //full x, every node
// 					for (int j = 0; j < nCAP; j++) //should be z
// 						MPI_Get(&slice[(i * nCAP + j) * n[0]], n[0], MPI_CXX_DOUBLE_COMPLEX, rank,
// 								(i * n[0] + (j + boxIndex * nCAP)) * n[0], n[0], MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);

// 								// from [0,sizeX]*n*n + ([0,nCAP]+box*nCAP)*n
// 			if (isTopX) //top
// 				for (int i = startX; i < n_lx[0]; i++)
// 					for (int j = 0; j < nCAP; j++)
// 						MPI_Get(&slice[(i * nCAP + j) * n[0]], n[0], MPI_CXX_DOUBLE_COMPLEX, rank,
// 								(i * n[0] + (j + boxIndex * nCAP)) * n[0], n[0], MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);

// 		}
// 		// n_lx[0],n_lx[1],nCAP
// 	}
// 	inline void getYBox(int rank, int boxIndex)
// 	{
// 		if constexpr (DIM == 2)
// 		{
// 			for (int i = 0; i < n_lx[0]; i++)
// 			{
// 				MPI_Get(&slice[i * n[0]], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
// 						i * n[0], nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 				MPI_Get(&slice[(i + 1) * n[0] - nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
// 						(i + 1) * n[0] - nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 			}

// 		}
// 		if constexpr (DIM == 3)
// 		{
// 			for (int i = 0; i < n_lx[0]; i++) //left
// 				for (int j = 0; j < nCAP; j++)
// 					MPI_Get(&slice[(i * n[0] + j) * nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
// 							(i * n[0] + j) * n[0] + boxIndex * nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 					//from [0, ]

// 			for (int i = 0; i < n_lx[0]; i++) //right
// 				for (int j = n[0] - nCAP; j < n[0]; j++)
// 					MPI_Get(&slice[(i * n[0] + j) * nCAP], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
// 							(i * n[0] + j) * n[0] + boxIndex * nCAP, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 						   //disp, count
// 		}
// 	}
// 	inline void getZBox(int rank, int boxIndex)
// 	{
// 		ind z_start;
// 		for (int i = 0; i < n_lx[0]; i++)
// 			for (int j = 0; j < nCAP; j++)
// 			{
// 				//front
// 				z_start = 0;
// 				MPI_Get(&slice[(i * nCAP + j) * n[0] + z_start], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
// 						(i * n[0] + (boxIndex * nCAP + j)) * n[0] + z_start, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 				//back
// 				z_start = n[0] - nCAP;
// 				MPI_Get(&slice[(i * nCAP + j) * n[0] + z_start], nCAP, MPI_CXX_DOUBLE_COMPLEX, rank,
// 						(i * n[0] + (boxIndex * nCAP + j)) * n[0] + z_start, nCAP, MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
// 				// disp, count
// 			}
// 	}
