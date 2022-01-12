#pragma once
#include <vector>
#include <mpi.h>
#include "fftw3.h"
#include "fftw3-mpi.h"
#include <map>

namespace MPI
{
	/* Global MPI communicator */
	int pID;			// Process id
	int pSize;			// # of processes
	MPI_Comm pComm = MPI_COMM_WORLD;
	MPI_Group pGroup;
	MPI_Status status;


	/* Current */
	int group;
	int region;
	int groupCount;
	int regionCount;
	/* Process id in MPI execution */
	int gID;	// In group
	int rID;	// In region
	int eID;	// In equiv. class
	/* Total number of MPI processes */
	int gSize;	// In group
	int rSize;	// In region
	int eSize;	// In equiv. class
	/* Communicate between */
	MPI_Comm gComm;
	MPI_Comm rComm;
	MPI_Comm eComm;
	MPI_Group gGroup;
	MPI_Group rGroup;
	MPI_Group eGroup;

	template <class T>
	constexpr auto type()
	{
		if constexpr (std::is_same_v<T, int>) return MPI_INT;
		if constexpr (std::is_same_v<T, double>) return MPI_DOUBLE;
		if constexpr (std::is_same_v<T, cxd>) return MPI_CXX_DOUBLE_COMPLEX;
		if constexpr (std::is_same_v<T, bool>) return MPI_CXX_BOOL;
	}

	void Barrier(MPI_Comm com = MPI_COMM_WORLD)
	{
		MPI_Barrier(com);
	}
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
	struct Division {};
	struct Slices :Division {};
	struct Rods :Division {};



	template <class MPIDivision, DIMS DIM, AXIS ...>
	struct Regions;

	template <DIMS DIM, AXIS ... freeAxes>
	struct Regions<Slices, DIM, freeAxes...>
	{
		using MPIDivision = Slices;
		static_assert(Power(2, DIM) >= sizeof...(freeAxes), "Number of regions must be less than 2^DIM");

		// Whether the current wf interacts in all spatial directions
		bool isMain;
		// Marks whether a direction is bounded (not free)
		bool bounded[DIM];
		AXIS freeCoord;
		int boundedCoordDim;

		static constexpr inline bool many = (sizeof...(freeAxes) > 1);
		static constexpr inline AXIS freeCoords[sizeof...(freeAxes)]
		{ freeAxes... };

		static constexpr inline int freeCoordsCount[sizeof...(freeAxes)]
		{ freeAxisCount<freeAxes>... };

		static inline int gMembers[uniq<freeAxisCount<freeAxes>...>::size];
		static inline int maxFree = 0;
		static inline int minFree = 100;

		static inline int groupLeader[uniq<freeAxisCount<freeAxes>...>::size];

		static inline MPI_Comm lessFree = 0; //region inter-comms
		static inline MPI_Comm moreFree = 0; //region inter-comms

		Regions()
		{
			logSETUP("Attempting to init %d groups and %d regions", groupCount, regionCount);
			regionCount = sizeof...(freeAxes);
			groupCount = uniq<freeAxisCount<freeAxes>...>::size;

			// assertm(MPI::rSize > regionCount, "Number of MPI processes too small");
			// logTestFatal(pSize % regionCount == 0, "Number of MPI processes pSize (%d) should be divisible by the number of regions regionCount (%d)", pSize, regionCount);

			//The following should work, but needs to be tested
			rSize = pSize / regionCount;
			region = pID / rSize; // Region number
			rID = pID - region * rSize;
			group = freeCoordsCount[region];
			freeCoord = freeCoords[region];
			boundedCoordDim = DIM - freeCoordsCount[region];
			isMain = (boundedCoordDim == DIM);
			for (DIMS i = 0; i < DIM; i++) //Determines whether individual directions are bounded (for ease of use)
				bounded[i] = !bool(freeCoord & getAxis(i));

			MPI_Comm_split(MPI_COMM_WORLD, group, pID, &gComm);
			MPI_Comm_rank(gComm, &gID);
			MPI_Comm_size(gComm, &gSize);
			MPI_Comm_group(gComm, &gGroup);
			_logMPI("[pID %3d/%3d] is in group %d [gID %3d/%3d]", pID, pSize, group, gID, gSize);

			Barrier();
			logSETUP("Initiated %d groups", groupCount);

			/* Other ideas: Interesting node local communicator:
			https://stackoverflow.com/questions/39912588/can-i-use-mpi-with-shared-memory */
			MPI_Comm_split(gComm, region, pID, &rComm);
			MPI_Comm_rank(rComm, &rID);
			MPI_Comm_size(rComm, &rSize);
			MPI_Comm_group(rComm, &rGroup);
			_logMPI("[pID %3d] region %d: [rID %3d/%3d] is in group %d [gID %3d/%3d]", pID, region, rID, rSize, group, gID, gSize);

			logSETUP("Initiated %d regions, processes per region: %d", regionCount, rSize);
			Barrier();

			//equivalence comm - for joining the wf
			MPI_Comm_split(MPI_COMM_WORLD, rID, pID, &eComm);
			MPI_Comm_rank(eComm, &eID);
			MPI_Comm_size(eComm, &eSize);
			MPI_Comm_group(eComm, &eGroup);

			_logMPI("pID %d belongs to equiv. class %d with eID %d/ %d", pID, rID, eID, eSize);
			Barrier();

			logSETUP("Attempting to init %d inter-communicators linking groups", groupCount - 1);
			// Determine number of members in each group (group is determined by freeCoordsCount)
			for (int i = 0; i < regionCount; i++)
				gMembers[freeCoordsCount[i]]++;

			// Determine (set) group p-leaders
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
			MPI_Comm lessFreeInter = 0; //region inter-comms
			MPI_Comm moreFreeInter = 0; //region inter-comms
			//inter-communicate with lessFree groups
			if (freeCoordsCount[region] > minFree)
			{
				MPI_Intercomm_create(gComm, 0, MPI_COMM_WORLD,
									 groupLeader[group - 1],
									 freeCoordsCount[region],
									 &lessFreeInter);

				MPI_Comm_remote_size(lessFreeInter, &size);
				// __logMPI("Process %d, local group size %d, Remote group (lessFree) size %d and should %d, minFree %d maxFree %d\n", pID, gMembers[freeCoordsCount[region]], size, gMembers[freeCoordsCount[region] - 1], minFree, maxFree);

				MPI_Intercomm_merge(lessFreeInter, 1, &lessFree);
			}
			//inter-communicate with moreFree groups
			if (freeCoordsCount[region] < maxFree)
			{
				MPI_Intercomm_create(gComm, 0, MPI_COMM_WORLD,
									 groupLeader[group + 1],
									 freeCoordsCount[region] + 1,
									 &moreFreeInter);

				MPI_Comm_remote_size(moreFreeInter, &size);
				// __logMPI("Process %d, local group size %d, Remote group (moreFree) size %d and should %d, minFree %d maxFree %d\n", pID, gMembers[freeCoordsCount[region]], size, gMembers[freeCoordsCount[region] + 1], minFree, maxFree);

				MPI_Intercomm_merge(moreFreeInter, 0, &moreFree);

				int rank;
				MPI_Comm_rank(moreFree, &rank);
				__logMPI("[group %d region %d] rank of process %d is %d\n", group, region, pID, rank);
			}
		}

		static bool calcsEnabled() { return group == 0; }
		static bool evoEnabled() { return true; };
	};

	template <DIMS D, class MPIDivision>
	using SingleRegion = MPI::Regions<MPIDivision, D, AXIS::NO>;

	// Here we merge positive and negative sides of axes into one region
	template <DIMS D, class MPIDivision>
	struct MultiRegionsReduced;

	template <class MPIDivision>
	struct MultiRegionsReduced<3_D, MPIDivision> : MPI::Regions<MPIDivision, 3_D,
		AXIS::NO, AXIS::X,
		AXIS::Y, AXIS::Z,
		AXIS::XY, AXIS::XZ,
		AXIS::YZ, AXIS::XYZ> {};

	template <class MPIDivision>
	struct MultiRegionsReduced<2_D, MPIDivision> : MPI::Regions<MPIDivision, 2_D,
		AXIS::NO, AXIS::X,
		AXIS::Y, AXIS::XY> {};

	template <class MPIDivision>
	struct MultiRegionsReduced<1_D, MPIDivision> : MPI::Regions<MPIDivision, 1_D,
		AXIS::NO, AXIS::X> {};


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

		// MPI_Scatter(mask, m_l,
		// 			is_same<T, double>() ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX,
		// 			local_v, m_l,
		// 			is_same<T, double>() ? MPI_DOUBLE : MPI_CXX_DOUBLE_COMPLEX,
		// 			0, MPI::rComm);
		// if (!MPI::rID) delete[] mask;
	}
};
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