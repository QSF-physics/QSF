#include "flux_mask.h"

borVec<DIM>* regionBorders;
template <typename Regions, typename FLUX_TYPE>
void defineBorder(size_t regionIndex, FLUX_TYPE flx)
{
	if (!MPI::pID)
	{
		for (i = 0; i < n; i++)
		{
			if constexpr (DIM == 1) mask[i] = flx.border(i - n2);
			else
			{
				readInd0 = i * n;
				for (j = 0; j < n; j++)
				{
					readInd1 = (readInd0 + j);
					if constexpr (DIM == 2)
						mask[readInd1] = flx.border(i - n2, j - n2);
					else
					{
						readInd1 = (readInd1)*n;
						for (k = 0; k < n; k++)
							mask[readInd1 + k] = flx.border(i - n2, j - n2, k - n2);
					}
				}
			}
		}
	#if EXTRA_OPTIONS & PREVIEW_REGIONS
		previewMask(mask, n, FLUX_TYPE::name);
	#endif
	}
	MPI_Barrier(MPI::rComm);
	MPI_Scatter(mask, local_m, borVec_MPI, &regionBorders[regionIndex * local_m],
				local_m, borVec_MPI, 0, MPI::rComm);
			// #if EXTRA_OPTIONS & PREVIEW_REGIONS_PER_NODE
			// 	logMPI("Border as distributed over nodes:");
			// 	for (int i = 0; i < MPI::rSize; i++)
			// 	{
			// 		if (MPI::pID == i)
			// 		{
			// 			previewMask(regionBorders[regionIndex], local_n, FLUX_TYPE::name);
			// 		}
			// 		MPI_Barrier(MPI::rComm);
			// 	}
			// #endif
}

template <typename PROP, typename... Args>
struct CO_FLUXES : COMPUTATION<double, Args...>, _RBUFFER
{
	static constexpr string_view name = "FLUX_";
	template <REP R, typename WHEN>
	constexpr static inline bool canRun = COMPUTATION<double, Args...>::template goodRep<R>;

	template <MODE M, REP R, OPTIMS opts> inline static void prepare()
	{
		logSETUP("Initializing flux %td regions [rowSize: %td]", (sizeof...(Args)), rowSize);
		p_NS = (ind)(12.5 / dx); 	// p_NS: position of borders between N and S regions (standard: 12.5)
		p_SD = (ind)(7. / dx);		// p_SD: position of borders between S and D regions (standard: 7)
		p_DT = (ind)(5. / dx);		// p_DT: position of borders between D and T regions (standard: 5)
		// p_CAP = (ind)(n2 - 1.0 * PROP::CAPlength / dx);// CAP_nodes;
		p_CAP = (ind)(20.0 / dx);// CAP_nodes;
		// logWarning("Number of cap nodes %g %td", PROP::CAPlength, ind(PROP::CAPlength / dx));
		logTest(p_NS > 0 && p_SD > 0, "Border radiuses p_NS (%td), p_SD (%td) positive", p_NS, p_SD);
		logTest(p_NS < p_CAP, "p_NS (%td) smaller than p_CAP (%td)", p_NS, p_CAP);
		curr = (borVec<DIM>*)malloc(sizeof(borVec<DIM>) * local_m);
		row_before = (cxd*)malloc(sizeof(cxd) * rowSize);
		row_after = (cxd*)malloc(sizeof(cxd) * rowSize);
		mask = (borVec<DIM>*)malloc(sizeof(borVec<DIM>) * m);
		regionBorders = new borVec<DIM>[(sizeof...(Args) * local_m)];
		// initMPIborVec(up_to<DIM - 1>);
		// MPI_Datatype borVec_MPI;
		MPI_Type_contiguous(DIM, MPI_DOUBLE, &borVec_MPI);
		MPI_Type_commit(&borVec_MPI);
		(defineBorder<typename PROP::Regions>(Index_v<Args, Args...>, Args{}), ...);
		MPI_Barrier(MPI::rComm);
		free(mask);
	}

	template <MODE M, REP R, OPTIMS opts> inline static void forerunner()
	{
		calcCurrentMap<PROP>();
	}
	template<REP R, OPTIMS opt, typename FLUX_TYPE> inline double calc() const noexcept
	{
		constexpr auto indx = Index_v<FLUX_TYPE, Args...>;
		size_t offset = indx * local_m;
		double res = 0.0;
		for (int i = 0; i < local_m; i++)
			res += (regionBorders[offset + i] * curr[i]);
		return res * dV / dx;
	}
};

// using ZOA_FLUX_3D = FLUX<N2S, N2D, N2T, S2D, S2T, S2CAP, D2CAP, T2CAP>;
// using ZOA_FLUX_2D = FLUX<N2S, N2D, S2D, S2CAP, D2CAP>;







//Using overload bordVec * bordVec = double and Green's theorem

// template <size_t ... args>
// void initMPIborVec(seq<args...> s)
// {
// 	int lengths[sizeof ... (args)] = { (0 * args + 1)... };
// 	const MPI_Aint displacements[sizeof ... (args)] = { (args * sizeof(double))... };
// 	MPI_Datatype types[sizeof ... (args)] = { (args >= 0 ? MPI_DOUBLE : MPI_DOUBLE)... }; //hackish
// 	MPI_Type_create_struct(sizeof ... (args), lengths, displacements, types, &borVec_MPI);
// 	MPI_Type_commit(&borVec_MPI);
// }

// void saveFluxRegions()
// {
	// #if RE_DATA & O_FLUXES
	// if (!MPI::pID)
	// {
	// 	file_flx = writeData(flx_prefix, "");
	// 	FILE* file_reg = writeData(reg_prefix, "");
	// 	double gs_pop = 1.0;
	// 	double single_pop = 0.0;
	// 	double double_pop = 0.0;
	// 	double absorbed_pop = 0.0;
	// 	for (int i = 0; i < scaled_ntsteps; i++)
	// 	{
	// 		fprintf(file_flx, "%i\t%.18E\t%.18E\t%.18E\t%.18E\t%.18E\t%.18E\n",
	// 				i * flux_interval, TOTAL_N2S[i], TOTAL_N2D_corr[i], TOTAL_N2D_anticorr[i], TOTAL_S2D[i], TOTAL_S2CAP[i], TOTAL_D2CAP[i]);

	// 		gs_pop -= (TOTAL_N2S[i] + TOTAL_N2D_corr[i] + TOTAL_N2D_anticorr[i]) * p_DT * flux_interval;
	// 		single_pop += (TOTAL_N2S[i] - TOTAL_S2D[i]) * p_DT * flux_interval;
	// 		double_pop += (TOTAL_N2D_corr[i] + TOTAL_N2D_anticorr[i]) * p_DT * flux_interval;
	// 		absorbed_pop += (TOTAL_S2CAP[i] + TOTAL_D2CAP[i]) * p_DT * flux_interval;
	// 		fprintf(file_reg, "%i\t%.18E\t%.18E\t%.18E\t%.18E\n",
	// 				i * flux_interval, gs_pop, single_pop, double_pop, absorbed_pop);
	// 	}

	// 	fclose(file_flx);
	// 	fclose(file_reg);
	// }
	// #endif
// }

