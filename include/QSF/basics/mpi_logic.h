#pragma once
#include <vector>
#include <mpi.h>
#include "fftw3.h"
#include "fftw3-mpi.h"
#include <map>

namespace MPI
{
	/* Current */
	int group;
	int region;

	int groupCount;
	int regionCount;
	MPI_Status status;
	/* Process id in MPI execution */
	int pID;	// Globally
	int gID;	// In group
	int rID;	// In region
	/* Total number of MPI processes */
	int pSize;	// Globally
	int gSize;	// In group
	int rSize;	// In region

	/* Communicate between */
	// MPI_COMM_WORLD;
	MPI_Comm gComm;
	MPI_Comm rComm;
	/* Groups of */
	MPI_Group pGroup;
	MPI_Group gGroup;
	MPI_Group rGroup;


	// All rigor opt: FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, FFTW_EXHAUSTIVE, FFTW_WISDOM_ONLY
	int plan_rigor = FFTW_MEASURE;

	char verstring[MPI_MAX_LIBRARY_VERSION_STRING];
	char nodename[MPI_MAX_PROCESSOR_NAME];
	int version, subversion, verstringlen, nodestringlen;
	void init(int argc, char* argv[])
	{
		MPI_Init(&argc, &argv);		 			 // initialize MPI
		MPI_Comm_rank(MPI_COMM_WORLD, &pID); 	// get id of this process
		MPI_Comm_size(MPI_COMM_WORLD, &pSize); // get number of processes
		MPI_Comm_group(MPI_COMM_WORLD, &pGroup);

		MPI_Get_processor_name(nodename, &nodestringlen);
		MPI_Get_version(&version, &subversion);
		MPI_Get_library_version(verstring, &verstringlen);

		if (pID == 0)
		{
			printf("MPI %d.%d support found\n", version, subversion);
			printf("MPI implementation: <%s>\n", verstring);
		}
	}

	void reduceImmediataly(double* variable, int size = 1)
	{
		MPI_Allreduce(MPI_IN_PLACE, variable, size, MPI_DOUBLE, MPI_SUM, rComm);
	}
	void reduceImmediataly(cxd* variable, int size = 1)
	{
		MPI_Allreduce(MPI_IN_PLACE, variable, size, MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, rComm);
	}

}


#define MPIB MPI_Barrier(MPI_COMM_WORLD);
#define _logMPI(...) { if constexpr (bool(DEBUG & DEBUG_MPI)) { for(int no=0; no<MPI::pSize; no++) { MPIB if(no==MPI::pID) _logColor(__LOG_YELLOW, __VA_ARGS__)  } }}
#define __logMPI(...) { if constexpr (bool(DEBUG & DEBUG_MPI)) { __LOG2S(__VA_ARGS__) }}

namespace MPI
{
	struct Strategy {};
	struct Slices :Strategy {};
	struct Rods :Strategy {};

	struct Single {};
	struct Multi {};

	struct SrcRegion
	{
		int rank;
		FREE_COORD fc;
	};

	std::vector<SrcRegion> sourceRegions;

	template <typename, DIMS D, FREE_COORD ...>
	struct Regions;

	template <DIMS D, FREE_COORD ... free_coord>
	struct Regions<Slices, D, free_coord...>
	{
		// static_assert(Power(2, DIM) >= sizeof...(free_coord), "Number of regions must be less than 2^DIM");
		using MPIStrategy = Slices;
		//Whether the current wf interacts in all spatial directions
		bool isMain;
		bool bounded[intDIMS(D)];
		FREE_COORD freeCoord;
		int boundedCoordDim;
		static constexpr inline bool many = (sizeof...(free_coord) > 1);
		static constexpr inline FREE_COORD freeCoords[sizeof...(free_coord)]
		{ free_coord... };

		static constexpr inline int freeCoordsCount[sizeof...(free_coord)]
		{ freeCoordCount<free_coord>... };

		static inline int gMembers[uniq<freeCoordCount<free_coord>...>::size];
		static inline int maxFree = 0;
		static inline int minFree = 100;

		static inline int groupLeader[uniq<freeCoordCount<free_coord>...>::size];


		static inline MPI_Comm lessFree = nullptr; //region inter-comms
		static inline MPI_Comm moreFree = nullptr; //region inter-comms

		Regions()
		{
			regionCount = sizeof...(free_coord);
			groupCount = uniq<freeCoordCount<free_coord>...>::size;
			logSETUP("Attempting to init %d groups and %d regions", groupCount, regionCount);
			// logTestFatal(pSize % regionCount == 0, "Number of MPI processes pSize (%d) should be divisible by the number of regions regionCount (%d)", pSize, regionCount);

			//The following should work, but needs to be tested
			rSize = pSize / regionCount;
			region = pID / rSize; // Region number
			rID = pID - region * rSize;
			group = freeCoordsCount[region];
			freeCoord = freeCoords[region];
			boundedCoordDim = intDIMS(D) - freeCoordsCount[region];
			// int boundedCoord = int(AXIS::ALL) - int(freeCoord);
			int boundedCoord = int(maxFreeCoord<D>) - int(freeCoord);
			isMain = (boundedCoordDim == intDIMS(D));
			for (int i = 0; i < intDIMS(D); i++)
			{
				//Determines whether X,Y,Z,... coords are bounded (for ease of use)
				bounded[i] = bool(boundedCoord & (1 << i));
			}

			MPI_Comm_split(MPI_COMM_WORLD, group, pID, &gComm);
			MPI_Comm_rank(gComm, &gID);
			MPI_Comm_size(gComm, &gSize);
			MPI_Comm_group(gComm, &gGroup);
			_logMPI("[pID %3d/%3d] is in group %d [gID %3d/%3d]", pID, pSize, group, gID, gSize);

			MPI_Barrier(MPI_COMM_WORLD);
			logSETUP("Initiated %d groups", groupCount);

			/* Other ideas: Interesting node local communicator:
			https://stackoverflow.com/questions/39912588/can-i-use-mpi-with-shared-memory */
			MPI_Comm_split(gComm, region, pID, &rComm);
			MPI_Comm_rank(rComm, &rID);
			MPI_Comm_size(rComm, &rSize);
			MPI_Comm_group(rComm, &rGroup);
			_logMPI("region %d: [rID %3d/%3d] is in group %d [gID %3d/%3d]", region, rID, rSize, group, gID, gSize);

			MPI_Barrier(MPI_COMM_WORLD);
			logSETUP("Initiated %d regions, processes per region: %d", regionCount, rSize);

			logSETUP("Attempting to init %d inter-communicators linking groups", groupCount - 1);
			// Determine number of members in each group
			for (int i = 0; i < regionCount; i++)
			{
				gMembers[freeCoordsCount[i]]++;
			}
			// Determine group p-leaders
			groupLeader[0] = 0;
			// logSETUP("Group %d has %d regions and p-leader %d", 0, gMembers[i], 0);
			for (int i = 1; i < groupCount; i++)
			{
				groupLeader[i] = groupLeader[i - 1] + rSize * gMembers[i - 1];
				logSETUP("Group %d has %d regions and p-leader %d (prev p-members %d)", i, gMembers[i], groupLeader[i], rSize * gMembers[i - 1]);
			}
			// Determine min/max group
			for (int i = 0; i < regionCount; i++)
			{
				if (freeCoordsCount[i] > maxFree)
					maxFree = freeCoordsCount[i];
				if (freeCoordsCount[i] < minFree)
					minFree = freeCoordsCount[i];
			}

			int size;
			MPI_Comm lessFreeInter = nullptr; //region inter-comms
			MPI_Comm moreFreeInter = nullptr; //region inter-comms
			//inter-communicate with lessFree groups
			if (freeCoordsCount[region] > minFree)
			{
				MPI_Intercomm_create(gComm, 0, MPI_COMM_WORLD,
									 groupLeader[group - 1], freeCoordsCount[region], &lessFreeInter);

				MPI_Comm_remote_size(lessFreeInter, &size);
				// __logMPI("Process %d, local group size %d, Remote group (lessFree) size %d and should %d\n", pID, gMembers[freeCoordsCount[region]], size, gMembers[freeCoordsCount[region] - 1]);

				MPI_Intercomm_merge(lessFreeInter, 1, &lessFree);
			}
			//inter-communicate with moreFree groups
			if (freeCoordsCount[region] < maxFree)
			{
				MPI_Intercomm_create(gComm, 0, MPI_COMM_WORLD,
									 groupLeader[group + 1], freeCoordsCount[region] + 1, &moreFreeInter);

				MPI_Comm_remote_size(moreFreeInter, &size);
				// __logMPI("Process %d, local group size %d, Remote group (moreFree) size %d and should %d\n", pID, gMembers[freeCoordsCount[region]], size, gMembers[freeCoordsCount[region] + 1]);

				MPI_Intercomm_merge(moreFreeInter, 0, &moreFree);
				int rank;
				MPI_Comm_rank(moreFree, &rank);
				// __logMPI("[group %d region %d] rank of process %d is %d\n", group, region, pID, rank);
			}

			int index = 0;
			for (int i = 0; i < regionCount; i++)
			{
				if (freeCoordsCount[region] - freeCoordsCount[i] == 1)
				{
					switch (int(freeCoord) - int(freeCoords[i]))
					{
					case (int)FREE_COORD::X: sourceRegions.push_back({ index * rSize + rID,FREE_COORD::X }); break;
					case (int)FREE_COORD::Y: sourceRegions.push_back({ index * rSize + rID,FREE_COORD::Y }); break;
					case (int)FREE_COORD::Z: sourceRegions.push_back({ index * rSize + rID,FREE_COORD::Z }); break;
					default: break;
					}
					index++;
				}
			}

			_logMPI("[group %d region %d pID %d] has %td sources", group, region, pID, sourceRegions.size());
		}

		static bool calcsEnabled()
		{
			return group == 0;
		}
		static bool evoEnabled()
		{
			return true;//region == 4;
		}
	};
}

template <DIMS D, class MPIStrategy>
struct MultiMPIGrid;

template <class MPIStrategy>
struct MultiMPIGrid<DIMS::D3, MPIStrategy> : MPI::Regions<MPIStrategy, DIMS::D3,
	FREE_COORD::NO, FREE_COORD::X, FREE_COORD::Y, FREE_COORD::Z, FREE_COORD::XY, FREE_COORD::XZ, FREE_COORD::YZ, FREE_COORD::XYZ>
{};

template <class MPIStrategy>
struct MultiMPIGrid<DIMS::D2, MPIStrategy> : MPI::Regions<MPIStrategy, DIMS::D2, FREE_COORD::NO, FREE_COORD::X, FREE_COORD::Y, FREE_COORD::XY> {};

template <class MPIStrategy>
struct MultiMPIGrid<DIMS::D1, MPIStrategy> :MPI::Regions<MPIStrategy, DIMS::D1, FREE_COORD::NO, FREE_COORD::X> {};

template <DIMS D, class MPIStrategy>
using SingleMPIGrid = MPI::Regions<MPIStrategy, D, FREE_COORD::NO>;

template <typename M, DIMS D, class MPIStrategy>
using MPIGrid = std::conditional_t < std::is_same_v<M, MPI::Single>, SingleMPIGrid<D, MPIStrategy>, MultiMPIGrid<D, MPIStrategy>>;

template <typename T, typename F>
void buildAndScatter(F& fun, T*& local_v)
{
	// T* mask;
	// if (!MPI::rID)
	// {
	// 	mask = new T[m];
	// 	fun(mask);
	// }
	// if (local_v) delete[] local_v;
	// local_v = new T[m];

	// MPI_Scatter(mask, local_m,
	// 			is_same<T, double>() ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX,
	// 			local_v, local_m,
	// 			is_same<T, double>() ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX,
	// 			0, MPI::rComm);
	// if (!MPI::rID) delete[] mask;
}

//MPI SIZE DEBUG
// ind loc0;
// 	ind loc1;
// 	ind locs0;
// 	ind locs1;
// 	fftw_mpi_local_size_2d_transposed(
// 		n, n, MPI_COMM_WORLD,
// 		&loc0, &locs0,
// 		&loc1, &locs1);
// 	logInfo("%td,%td,%td,%td", loc0, locs0, loc1, locs1);