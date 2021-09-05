namespace QSF
{
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
		using Base::m_l, Base::m, Base::xmin, Base::dx, Base::strides_lx;
		using Base::mcomm, Base::psi, Base::mask, Base::inv_m, Base::inv_n;
		static_assert(Base::hasAbsorber, "Absorber is required for multigrid computations");
		static inline constexpr cxd zero = { 0.0, 0.0 };
		template <ind Is> ind static constexpr rev = DIM - 2 - Is; //when iterating over n-1 dims

		MPI_Win lessFreeWin = NULL;
		MPI_Win moreFreeWin = NULL;
		bool mpiFFTW;
		fftw_plan extra_plans[2 * DIM];
		int extra_plans_count = 0;

		ind stride_slice[DIM][DIM]; //For each "free" axis
		fftw_plan transf_slice[DIM];
		cxd* slice;
		ind m_slice;
		ind m_lslice;
		int boxesCount;
		bool isBottomX;
		bool isTopX;
		ind sizeX;
		ind startX;

		std::vector<SrcRegion> sourceRegions;

		LocalGrid(Base g) : Base(g), mpiFFTW(g.mcomm.bounded[0]) { init(); }
		LocalGrid(Section& settings) :Base(settings), mpiFFTW(Base::mcomm.bounded[0]) { init(); }

		void initExtraFFTW()
		{
			if (!mcomm.isMain)
			{

				if (mpiFFTW)
				{
					printf("pID %d: is not main makes mpiFFTW plans\n", MPI::pID);
					//make forward and backward mpi plans
					for (int i = 0; i < 2; i++)
						Base::mpi_plans[i] =
						fftw_mpi_plan_many_dft(1, n, m / n[0], //rank, dims, howmany
											   FFTW_MPI_DEFAULT_BLOCK, //block 
											   FFTW_MPI_DEFAULT_BLOCK, //tblock
											   reinterpret_cast<fftw_complex*>(psi),
											   reinterpret_cast<fftw_complex*>(psi),
											   MPI::rComm, Base::for_back[i], MPI::plan_rigor);


				}

				int fftw_rank = 0;
				int start_dim = -1;
				extra_plans_count = 0;
				// Non-mpi fftw (extra) plans for non-main region (region!=0) are coupled 
				// together if involve consecutive directions (case for 1-3D)

				for (DIMS i = 1; i <= DIM; i++)
				{
					logWarning("FFTW extra steps");
					if (i == DIM || !mcomm.bounded[i])
					{
						if (fftw_rank)
						{
							fftw_iodim64* dims = new fftw_iodim64[fftw_rank];
							for (int j = 0; j < fftw_rank; j++)
							{
								dims[j] = { n_lx[start_dim + j], strides_lx[start_dim + j],strides_lx[start_dim + j] };
								// __logMPI("region %d start_dim %d j %d n_lx[start_dim+j]=%td strides_lx[start_dim+j]=%td\n", MPI::region, start_dim, j, n_lx[start_dim + j], strides_lx[start_dim + j]);
							}

							//define a loop over all howmany_dims (free) dimensions
							int index = 0;
							int howmany_rank = DIM - fftw_rank;
							fftw_iodim64* howmany_dims = new fftw_iodim64[howmany_rank];
							for (int j = 0; j < DIM;j++)
							{
								if (j < start_dim || (j >= start_dim + fftw_rank && j < DIM))
								{
									howmany_dims[index] = { n_lx[j], strides_lx[j], strides_lx[j] };
									index++;
									// __logMPI("region %d start_dim %d index %d j %d n_lx[j]=%td dist[j]=%td\n", MPI::region, start_dim, index, j, n_lx[j], strides_lx[j]);
								}
							};
							for (int k = 0; k < 2; k++)
								extra_plans[extra_plans_count + DIM * k] = fftw_plan_guru64_dft(
									fftw_rank, dims,
									howmany_rank, howmany_dims,
									reinterpret_cast<fftw_complex*>(psi),
									reinterpret_cast<fftw_complex*>(psi),
									Base::for_back[k], MPI::plan_rigor);

							delete[] dims;
							delete[] howmany_dims;
							extra_plans_count++;
							// __logMPI("region %d mpiFFTW %d extra_plans_count %d dim %d/%d fftw_rank %d howmany_rank %d\n", MPI::region, mpiFFTW, extra_plans_count, start_dim, i - 1, fftw_rank, howmany_rank);
						}
						fftw_rank = 0;
					}
					else if (mcomm.bounded[i])
					{
						if (start_dim == -1) start_dim = i;
						fftw_rank++;
					}
				}
				Base::reset();
			}
			// fprintf(stderr, "---> region %d mpiFFTW %d extra_plans_count %d\n", MPI::region, mpiFFTW, extra_plans_count);
			// _logMPI("Region %d MPI::pID: %d got %td nodes or %td rows [start row: %td end row: %td]", MPI::region, MPI::pID, m_l, n_lx[0], pos_lx.first, pos_lx.last);
		}

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
			if (!mcomm.isMain)
			{
				initExtraFFTW();
				initSlice(n_seq<DIM>); //Main region doesn't operate with slices

			}
			if (mcomm.moreFree) MPI_Win_create(psi, m_l * sizeof(cxd), sizeof(cxd), MPI_INFO_NULL, mcomm.moreFree, &moreFreeWin);
			if (mcomm.lessFree) MPI_Win_create(MPI_BOTTOM, 0, sizeof(cxd), MPI_INFO_NULL, mcomm.lessFree, &lessFreeWin);
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
				// printf("pID %d MPI initFFTWplans dirFree: [%td] howmany %td\n",
					//    MPI::pID, dirFree, m_slice / n[0]);
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
				// printf("pID %d initFFTWplans dirFree: [%td] howmany %td\n", MPI::pID, dirFree, m_slice / n[0]);
				// The dimension we're going to transform
				fftw_iodim64 dims[1] = { {n[dirFree] , stride_slice[dirFree][dirFree], stride_slice[dirFree][dirFree]} };
				fftw_iodim64 howmany_dims[DIM - 1];
				int index = 0;
				// if (MPI::pID == 1)
					// std::cout << "---:> pID:" << MPI::pID << " dirFree:" << dirFree << " dims: " << dims[0].n << " " << dims[0].is << std::endl;

				([&] {if (dirs != dirFree) {
					howmany_dims[index] = { sliceShape<dirFree, dirs>(), stride_slice[dirFree][dirs], stride_slice[dirFree][dirs] };
					index++;
					// if (MPI::pID == 1)
						// std::cout << "---> pID:" << MPI::pID << " dirFree:" << dirFree << " dirs:" << dirs << "howmany_dims: " << howmany_dims[index - 1].n << " " << howmany_dims[index - 1].is << std::endl;
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

			printf("pID: %d  startX[%td] (sizeX[%td] isBottomX %d) isTopX %d \n", MPI::pID, startX, sizeX, isBottomX, isTopX);
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
			if (slice_off<0 || slice_off + count >m_lslice)
				printf("ERROR %d slice[%td/%td/%td]", MPI::pID, slice_off, count, m_lslice);
			if (wf_offset<0 || wf_offset + count >m_l)
				printf("ERROR %d wf[% td / % td / % td] \n", MPI::pID, wf_offset, count, m_l);
			MPI_Get(&slice[slice_off], //pos in slice
					count, //send count
					MPI_CXX_DOUBLE_COMPLEX, rank,
					wf_offset, //pos in wf
					count, //recv count
					MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
		}
		void backup(const ind step)
		{
			Base::_backup(step, IO::psi_ext + std::to_string(MPI::region));
		}
		void restore()
		{
			Base::_load("backup", IO::psi_ext + std::to_string(MPI::region));
			logWarning("Restored the wf for all regions");
		}
		void remove_backups()
		{
			Base::remove_backups(IO::psi_ext + std::to_string(MPI::region));
		}
		std::string save(std::string common_name = "", DUMP_FORMAT df = { .dim = DIM })
		{
			if (df.rep == REP::P) fourier<REP::P>();
			return Base::_save(common_name, df, IO::psi_ext + std::to_string(MPI::region));
			if (df.rep == REP::P) fourier<REP::X>();
		}

		std::string saveJoined(std::string common_name = "", DUMP_FORMAT df = { .dim = DIM })
		{
			if (df.rep == REP::P) fourier<REP::P>();
			// MPI_Allreduce(MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI::eComm);
			MPI_Reduce((MPI::eID) ? psi : MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI::eComm);
			if (MPI::region == 0) return Base::_save(common_name + "_joined", df);
			if (df.rep == REP::P) fourier<REP::X>();
		}

		std::string saveIonizedJoined(std::string common_name = "", DUMP_FORMAT df = { .dim = DIM })
		{
			if (df.rep == REP::P) fourier<REP::P>();
			// MPI_Allreduce(MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI::eComm);
			if (MPI::region == 0) Base::reset(); //removing un-ionized part

			MPI_Reduce((MPI::eID) ? psi : MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI::eComm);
			return Base::_save(common_name + "_ionized_joined", df);
			if (df.rep == REP::P) fourier<REP::X>();
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

		// inline static constexpr bool corrections = true;

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
			if constexpr (DIM == 3)
				return  onehalf * ((dirs == dirFree ? 0 : Base::template mask<dirs>(nodes)) + ... + 0)
				+ onesixth * (1 - ((dirs == dirFree ? 1.0 : Base::template mask<dirs>(nodes))*...));
			// else if constexpr (DIM == 2) //TODO: FIX
				// return onehalf * (1 - ((dirs == dirFree ? 1.0 : Base::template mask<dirs>(nodes))*...));
			else return 1.0;
				// else return 1.0;
			// onesixth* (invCAP2(j + boxIndex * CAPnodes, k)) + onehalf * (CAP1(j + boxIndex * CAPnodes) + CAP1(k));
		}

		template <uind dirFree, uind ... dirs>
		inline void maskSlice(int boxIndex, seq<dirs...>)
		{
			ind counters[DIM] = { 0 };
			do
			{
				//those by design operate on nCAP nodes, which takes neg dist from boundry
				if (Axis<dirFree> != AXIS::X || (isBottomX && counters[0] < sizeX))
				{
					slice[slice_offset<dirFree>(counters[dirs]...)]
						*= Base::template inv_mask<0>(Base::template neg_dist_from_edge<dirFree>(counters[dirFree]));

					if (MPI::group == 1)
						slice[slice_offset<dirFree>(counters[dirs]...)]
						*= mask_correction<dirFree, dirs...>(Base::template neg_dist_from_edge<dirs>(counters[dirs] + slider<dirFree, dirs>(boxIndex))...);
					// if ((MPI::group == 1) && (dirFree == 0))// && ((counters[dirs] == 1) && ...))
					// 	printf("pID %d reg %d box %d dir %td cooords[%2td %2td %2td] local[%2td %2td %2td] (bottom) val [%g]\n", MPI::pID, MPI::group, boxIndex, dirFree, Base::template neg_dist_from_edge<dirs>(counters[dirs] + slider<dirFree, dirs>(boxIndex))...,
					// 		   counters[dirs]...,
					// 		   mask_correction<dirFree, dirs...>(Base::template neg_dist_from_edge<dirs>(counters[dirs] + slider<dirFree, dirs>(boxIndex))...));
				}

				if (Axis<dirFree> != AXIS::X || (isTopX && counters[0] >= startX))
				{
					slice[slice_offset<dirFree>((counters[dirs] + top<dirFree, dirs>())...)]
						*= Base::template inv_mask<0>(Base::template neg_dist_from_edge<dirFree>(counters[dirFree] + top<dirFree, dirFree>()));
					if (MPI::group == 1)
						slice[slice_offset<dirFree>((counters[dirs] + top<dirFree, dirs>())...)]
						*= mask_correction<dirFree, dirs...>(Base::template neg_dist_from_edge<dirs>(counters[dirs] + tslider<dirFree, dirs>(boxIndex))...);

					// if ((MPI::group == 1) && (dirFree == 0))// && ((counters[dirs] == 1) && ...))
					// 	printf("pID %d reg %d box %d dir %td cooords[%2td %2td %2td] local[%2td %2td %2td] (top   ) val [%g]\n", MPI::pID, MPI::group, boxIndex, dirFree, Base::template neg_dist_from_edge<dirs>(counters[dirs] + tslider<dirFree, dirs>(boxIndex))...,
					// 		   counters[dirs]...,
					// 		   mask_correction<dirFree, dirs...>(Base::template neg_dist_from_edge<dirs>(counters[dirs] + tslider<dirFree, dirs>(boxIndex))...));
				}
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
			if (MPI::group != MPI::groupCount - 1)
				Base::maskRegion();
		}
		const unsigned MPI_win_flags = MPI_MODE_NOSTORE | MPI_MODE_NOPUT;
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
						if (MPI::group == gI - 1) MPI_Win_fence(MPI_win_flags | MPI_MODE_NOPRECEDE, moreFreeWin);
						if (MPI::group == gI)
						{
							//TODO: check for null ref?
							auto& src = sourceRegions[sI];
							resetSlice();
							MPI_Win_fence(MPI_win_flags | MPI_MODE_NOPRECEDE, lessFreeWin);
							boxDispatcher(src.fc, src.rank, bI);
							MPI_Win_fence(MPI_win_flags | MPI_MODE_NOSUCCEED, lessFreeWin);
							coherentAddition(src.fc, bI);
						}
						if (MPI::group == gI - 1) MPI_Win_fence(MPI_win_flags | MPI_MODE_NOSUCCEED, moreFreeWin);
					}
				}
			}
			// Timings::measure::stop("TRANSFER");
		}

	#pragma region FFToverloads
		inline void normalizeAfterTwoFFT()
		{
			if constexpr (Base::MPIGridComm::many)
			{
				if (mcomm.boundedCoordDim > 0)
					Base::multiplyArray(psi, Power(inv_n[0], mcomm.boundedCoordDim));
			}
			else Base::normalizeAfterTwoFFT();
		}


		template <REP R>
		inline void fourier()
		{
			constexpr DIMS back = R == REP::P ? 0 : 1;

			if (mpiFFTW)
				fftw_execute(Base::mpi_plans[back]);

			for (int i = 0; i < extra_plans_count; i++)
				fftw_execute(extra_plans[i + DIM * back]);

			if constexpr (back)
				normalizeAfterTwoFFT();
		}

	#pragma endregion FFToverloads

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
};