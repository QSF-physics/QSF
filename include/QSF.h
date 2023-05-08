#include "basics/constants.h"
#include "basics/data_types.h"
#include "basics/macros.h"
#include <filesystem>
// #include "basics/globals.h" //delete!
#include "basics/config.h"
#include "basics/enums.h"
#include "basics/inipp.h"
#include "basics/io.h"
#include "basics/logging.h"
#include "basics/math_basic.h"
#include "basics/mpi_logic.h"
#include "basics/static_seq.h"
#include "basics/static_string.h"
#include "basics/static_utils.h"
#include "basics/superpos.h"
#include "basics/timings.h"

#include "fluxes/borders.h"
#include "wf/absorbers.h"
#include "wf/buffer.h"
#include "wf/computations.h"
#include "wf/coords.h"
#include "wf/grid.h"
#include "wf/multigrid.h"
#include "wf/preset.h"
// #include "dumps.h"
// #include "average.h"
#include "coupling.h"
#include "field.h"
#include "hamiltonian.h"
#include "potential.h"
#include "propagator.h"
#include "routines.h"
namespace QSF
{
	void
	init(int argc, char* argv[], std::filesystem::path location= IO::project_dir / IO::results_dir)
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

	void subdirectory(std::filesystem::path sub)
	{
		if(!MPI::pID) bool success= std::filesystem::create_directories(IO::root_dir / sub);

		MPI::Barrier();
		std::filesystem::current_path(IO::root_dir / sub);
	}
	void finalize() { MPI_Finalize(); }
}		// namespace QSF