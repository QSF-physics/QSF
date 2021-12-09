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
		//Base (non-multi) grid type
		using Base = LocalGrid<Hamiltonian, BaseGrid, Components, MPI_GC, MPI::Slices, false>;
		using Base::DIM, Base::n_lx, Base::n, Base::pos_lx, Base::nCAP;
		using Base::m_l, Base::m, Base::xmin, Base::dx, Base::strides_lx;
		using Base::mcomm, Base::psi, Base::mask, Base::inv_m, Base::inv_n;
		static_assert(Base::hasAbsorber, "Absorber is required for multigrid computations");

		//Reverse indices order. Used when iterating over n-1 dims
		template <ind Is> ind static constexpr rev = DIM - 2 - Is;

		MPI_Win lessFreeWin = NULL;
		MPI_Win moreFreeWin = NULL;
		bool mpiFFTW;
		fftw_plan mpi_plans_partial[2];
		fftw_plan plans_partial[2 * DIM];
		int plans_partial_count = 0;
		double partial_inv_m;

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

		void initFFTW_partial()
		{
			if (!mcomm.isMain)
			{
				if (mpiFFTW)
				{
					printf("PARTIAL_FFTW[0]: pID %d: makes mpiFFTW plans of rank 1, howmany %td\n", MPI::pID, m / n[0]);
					//make forward and backward mpi plans
					for (int i = 0; i < 2; i++)
						mpi_plans_partial[i] =
						fftw_mpi_plan_many_dft(1, n, m / n[0], //rank, dims, howmany
											   FFTW_MPI_DEFAULT_BLOCK, //block 
											   FFTW_MPI_DEFAULT_BLOCK, //tblock
											   reinterpret_cast<fftw_complex*>(psi),
											   reinterpret_cast<fftw_complex*>(psi),
											   MPI::rComm, Base::for_back[i], MPI::plan_rigor);
				}

				int fftw_rank = 0;
				int start_dim = -1;
				plans_partial_count = 0;
				// Non-mpi fftw (partial) plans for non-main region (region!=0) are coupled 
				// together if involve consecutive directions (case for 1-3D)
				for (DIMS i = 1; i <= DIM; i++)
				{
					if (i == DIM || !mcomm.bounded[i])
					{
						if (fftw_rank)
						{
							fftw_iodim64* dims = new fftw_iodim64[fftw_rank];
							for (int j = 0; j < fftw_rank; j++)
							{
								dims[j] = { n_lx[start_dim + j], strides_lx[start_dim + j],strides_lx[start_dim + j] };
								// __logMPI("PARTIAL_FFTW[1]: pID %d region %d start_dim %d j %d n_lx[start_dim+j]=%td strides_lx[start_dim+j]=%td\n", MPI::pID, MPI::region, start_dim, j, n_lx[start_dim + j], strides_lx[start_dim + j]);
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
									// __logMPI("PARTIAL_FFTW[2]: pID %d region %d start_dim %d index %d j %d n_lx[j]=%td dist[j]=%td\n", MPI::pID, MPI::region, start_dim, index, j, n_lx[j], strides_lx[j]);
								}
							};
							for (int k = 0; k < 2; k++)
								plans_partial[plans_partial_count + DIM * k] = fftw_plan_guru64_dft(
									fftw_rank, dims,
									howmany_rank, howmany_dims,
									reinterpret_cast<fftw_complex*>(psi),
									reinterpret_cast<fftw_complex*>(psi),
									Base::for_back[k], MPI::plan_rigor);

							delete[] dims;
							delete[] howmany_dims;
							plans_partial_count++;
							// __logMPI("PARTIAL_FFTW[3]: pID %d region %d mpiFFTW %d plans_partial_count %d dim %d/%d fftw_rank %d howmany_rank %d\n", MPI::pID, MPI::region, mpiFFTW, plans_partial_count, start_dim, i - 1, fftw_rank, howmany_rank);
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
			// fprintf(stderr, "---> region %d mpiFFTW %d plans_partial_count %d\n", MPI::region, mpiFFTW, plans_partial_count);
			// _logMPI("Region %d MPI::pID: %d got %td nodes or %td rows [start row: %td end row: %td]", MPI::region, MPI::pID, m_l, n_lx[0], pos_lx.first, pos_lx.last);
		}

		void init()
		{
			__logMPI("MultiRegions init pID:%d", MPI::pID);
			int index = 0;
			for (int i = 0; i < MPI::regionCount; i++)
			{
				if (mcomm.freeCoordsCount[MPI::region] - mcomm.freeCoordsCount[i] == 1)
				{
					switch (mcomm.freeCoord ^ mcomm.freeCoords[i])
					{
					case AXIS::X: sourceRegions.push_back({ index * MPI::rSize + MPI::rID,AXIS::X });
						// printf("[[[%d is X source for %d]]]\n", index * MPI::rSize + MPI::rID, MPI::pID);
						break;
					case AXIS::Y: sourceRegions.push_back({ index * MPI::rSize + MPI::rID,AXIS::Y });
						// printf("[[[%d is Y source for %d]]]\n", index * MPI::rSize + MPI::rID, MPI::pID);
						break;
					case AXIS::Z: sourceRegions.push_back({ index * MPI::rSize + MPI::rID,AXIS::Z });
						// printf("[[[%d is Z source for %d]]]\n", index * MPI::rSize + MPI::rID, MPI::pID);
						break;
					default: break;
					}
					index++;
				}
			}
			//Used for normalization after partial FFT done two times
			partial_inv_m = 1.0;
			for (DIMS i = 0; i < DIM; i++)
				partial_inv_m *= mcomm.bounded[i] ? inv_n[i] : 1.0;


			// _logMPI("[group %d region %d pID %d] has %td sources, [isMain:%d]", MPI::group, MPI::region, MPI::pID, sourceRegions.size(), mcomm.isMain);

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

			_logMPI("boxesCount %d nCAP %td", boxesCount, nCAP);



			initFFTW_partial();
			initSlice(n_seq<DIM>); //Main region doesn't operate with slices



			if (mcomm.lessFree) MPI_Win_create(MPI_BOTTOM, 0, sizeof(cxd), MPI_INFO_NULL, mcomm.lessFree, &lessFreeWin);

			if (mcomm.moreFree) MPI_Win_create(psi, m_l * sizeof(cxd), sizeof(cxd), MPI_INFO_NULL, mcomm.moreFree, &moreFreeWin);
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
			// std::cout << MPI::pID << " dirFree: " << dirFree << " stride_slice[dirFree][..]" << stride_slice[dirFree][0] << " " << stride_slice[dirFree][1] << " " << stride_slice[dirFree][2] << " " << std::endl;
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
		void initSliceFFTWplans(seq<dirs...>)
		{
			if constexpr (dirFree == 0)
			{
				printf("SLICE_FFTW[0] pID %d MPI initSliceFFTWplans dirFree: [%td] howmany %td\n",
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
				printf("SLICE_FFTW[1] pID %d initSliceFFTWplans dirFree: [%td] howmany %td\n", MPI::pID, dirFree, m_slice / n[0]);
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
			if (!mcomm.isMain)
			{
				m_slice = DIM == 3 ? m / n[2] * nCAP : m;
				m_lslice = m_slice / MPI::rSize;
				slice = (cxd*)fftw_malloc(sizeof(cxd) * m_lslice);
				// __logMPI("[pID:%d] initSlice [m_slice:%td] [m_lslice:%td]\n", MPI::pID, m_slice, m_lslice);

				// Base::initFFTW();
				(initSliceStrides<dirFree>(rev_seq<DIM - 1>), ...);
				(initSliceFFTWplans<dirFree>(n_seq<DIM>), ...);

				isBottomX = pos_lx.first < nCAP;
				sizeX = Min(n_lx[0], nCAP - pos_lx.first);

				isTopX = pos_lx.last > n[0] - 1 - nCAP;
				startX = Max(0, n[0] - nCAP - pos_lx.first);

				// __logMPI("[xSLICE] pID: %d  startX[%td] (sizeX[%td] isBottomX %d) isTopX %d \n", MPI::pID, startX, sizeX, isBottomX, isTopX);
				if (transf_slice[0] == NULL) logInfo("fftw_mpi_plan... (slice) returned NULL");
			}
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
			else
				return sliceShape<dirFree, dir>();
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
			// static int a = 0;
			// if (a < 1000)
			// 	__logMPI("[pID:%d [%d]] [slice_off:%td] [count:%td] [rank:%d] [wf_offset:%td]\n", MPI::pID, a, slice_off, count, rank, wf_offset);
			// a++;
			MPI_Get(&slice[slice_off], //pos in slice
					count, //send count
					MPI_CXX_DOUBLE_COMPLEX, rank,
					wf_offset, //pos in wf
					count, //recv count
					MPI_CXX_DOUBLE_COMPLEX, lessFreeWin);
		}
		void backup(const ind step)
		{
			logIO("Multiregion backup");
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

		// On multiregion grid some directions are never transformed to X representation
		// hence in order to save them in homogenous (P or X) rep we first need to transform
		// them to P rep using multigrid::fourier_partial, then to X using (single) grid::fourier
		// and leave the rest of work to grid::snapshot
		void snapshot(std::string extra_info, DUMP_FORMAT df = { .dim = DIM })
		{
			logIO("multireg snapshot");
			if (!mcomm.isMain)
			{
				fourier<REP::P>();
				Base::template fourier<REP::X>();
			}
			Base::snapshot(extra_info, df);
			if (!mcomm.isMain)
			{
				Base::template fourier<REP::P>();
				fourier<REP::X>();
			}
		}

		std::string save(std::string common_name = "", DUMP_FORMAT df = { .dim = DIM })
		{
			if (!mcomm.isMain)
			{
				fourier<REP::P>();
				Base::template fourier<REP::X>();
			}
			auto ret = Base::_save_centered(common_name, df, IO::psi_ext + std::to_string(MPI::region));
			if (!mcomm.isMain)
			{
				Base::template fourier<REP::P>();
				fourier<REP::X>();
			}
			return ret;
		}

		// std::string saveJoined(std::string common_name = "", DUMP_FORMAT df = { .dim = DIM })
		// {
		// 	if (df.rep == REP::P) fourier<REP::P>();

		// 	// MPI_Allreduce(MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI::eComm);
		// 	MPI_Reduce((MPI::eID) ? psi : MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI::eComm);
		// 	if (MPI::region == 0) return Base::_save(common_name + "_joined", df);
		// 	if (df.rep == REP::P) fourier<REP::X>();
		// }

		std::string saveIonizedJoined(std::string common_name = "", DUMP_FORMAT df = { .dim = DIM })
		{
			if (!mcomm.isMain)
			{
				fourier<REP::P>();
				Base::template fourier<REP::X>();
			}
			if (MPI::region == 0) Base::reset(); //removing un-ionized part

			// MPI_Allreduce(MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI::eComm);
			MPI_Reduce((MPI::eID) ? psi : MPI_IN_PLACE, psi, m_l, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, 0, MPI::eComm);

			//EXPERIMENTAL
			IO::path ret = "";
			if (MPI::region == 0)
			{
				ret = Base::_save_centered(common_name + "_ionized_joined", df);
			}
			if (!mcomm.isMain)
			{
				Base::template fourier<REP::P>();
				fourier<REP::X>();
			}
			return ret;
		}


		template <uind dirFree, uind ... dirs>
		inline void getBox(int rank, int box, seq<dirs...>)
		{
			// return;
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
			switch (fc)
			{
			case AXIS::X:
				getBox<0>(rank, boxIndex, n_seq<DIM - 1>);
				break;
			case AXIS::Y:
				getBox<1>(rank, boxIndex, n_seq<DIM - 1>);
				break;
			case AXIS::Z:
				getBox<2>(rank, boxIndex, n_seq<DIM - 1>);
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
				// return  onehalf * ((dirs == dirFree ? 0 : Base::template mask<dirs>(nodes)) + ... + 0);
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
		inline void normalize_after_fft_partial()
		{
			if (mcomm.boundedCoordDim > 0)
				Base::multiplyArray(psi, partial_inv_m);
		}

		template <REP R>
		inline void fourier_partial()
		{
			// printf("multigrid FT %d pID:%d mpiFFTW:%d plans_partial_count:%d isMAIN:%d\n", int(R), MPI::pID, mpiFFTW, plans_partial_count, mcomm.isMain);
			constexpr DIMS back = R == REP::P ? 0 : 1;
			if (mpiFFTW) fftw_execute(mpi_plans_partial[back]);
			for (int i = 0; i < plans_partial_count; i++) fftw_execute(plans_partial[i + DIM * back]);
			if constexpr (back) normalize_after_fft_partial();
		}
		template <REP R>
		inline void fourier()
		{
			if (mcomm.isMain) Base::template fourier<R>();
			else fourier_partial<R>();
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