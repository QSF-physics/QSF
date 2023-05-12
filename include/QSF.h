/// NOTE - Do not change the order of includes!
#include <filesystem>
#include "basics/macros.h"
#include "basics/data_types.h"
#include "basics/constants.h"
// #include "basics/globals.h" //delete!
#include "basics/enums.h"
#include "basics/math_basic.h"
#include "basics/logging.h"
#include "basics/static_seq.h"
#include "basics/static_string.h"
#include "basics/static_utils.h"
#include "basics/mpi_logic.h"
#include "basics/superpos.h"
#include "basics/io.h"
#include "basics/timings.h"
#include "basics/inipp.h"
#include "basics/config.h"

#include "fluxes/borders.h"
#include "wf/preset.h"
#include "wf/coords.h"
#include "wf/absorbers.h"
#include "wf/computations.h"
#include "wf/buffer.h"
#include "wf/grid.h"
#include "wf/multigrid.h"
// #include "dumps.h"
// #include "average.h"
#include "field.h"
#include "coupling.h"
#include "hamiltonian.h"
#include "potential.h"
#include "propagator.h"
namespace QSF
{

	/// @brief Initializes the environment for a parallel program with MPI and sets up the project
	/// directory.
	/// @param argc The count of command line arguments
	/// forwarded to the MPI::init
	/// @param argv The vector of command line argument
	/// forwarded to the MPI::init
	/// @param location  The filesystem path for the root directory of the project
	// If not provided, the default is the results directory inside the project directory.
	void init(
		int argc, char* argv[], std::filesystem::path location= IO::project_dir / IO::results_dir)
	{
		MPI::init(argc, argv);
		IO::root_dir= location;
		if(!MPI::pID) bool success= std::filesystem::create_directories(location);
		MPI::Barrier();
		std::filesystem::current_path(location);
		// We forward argc, argv arguments, but this is NOT part of the MPI standard!
		// See W. Gropp et al. - Using MPI Portable Parallel Programming with the Message-Passing
		// Interface (2014, The MIT Press), p.60
		logImportant("PROJECT: [%s] MPI PROCESSES: [%d]", IO::project_name.c_str(), MPI::pSize);
		logImportant("MAIN OUTPUT PATH: [%s]", location.c_str());
	}
	/// @brief Calling this will redirect any futher program output to a subdirectory
	/// @param sub_path subdirectory tree
	void subdirectory(std::filesystem::path sub_path)
	{
		if(!MPI::pID) bool success= std::filesystem::create_directories(IO::root_dir / sub_path);

		MPI::Barrier();
		std::filesystem::current_path(IO::root_dir / sub_path);
	}
	/// @brief Called at the end of each cpp project file to gracefully
	/// exit the multiprocess program
	void finalize() { MPI_Finalize(); }
}		// namespace QSF