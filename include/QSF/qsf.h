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


#include "wf/absorbers.h"
#include "fluxes/borders.h"
#include "wf/preset.h"
#include "wf/coords.h"
#include "wf/computations.h"
#include "wf/buffer.h"
#include "wf/grid.h"
#include "wf.h"
// #include "dumps.h"
// #include "average.h"
#include "coupling.h"
#include "hamiltonian.h"
#include "potential.h"
#include "field.h"
#include "propagator.h"
#include "routines.h"

// #include "setup/buffer.h"
namespace QSF
{
	void init(int argc, char* argv[])
	{
		// We forward argc, argv arguments, but this is NOT part of the MPI standard!
		// See W. Gropp et al. - Using MPI Portable Parallel Programming with the Message-Passing Interface (2014, The MIT Press), p.60
		MPI::init(argc, argv);
		logImportant("PROJECT: [%s] MPI PROCESSES: [%d]", STRINGIFY(PROJNAME), MPI::pSize);
		createDir(project_dir, results_dir);
		// createDir(home_dir, results_project_dir);
	}

	void finalize()
	{
		MPI_Finalize();
	}
}