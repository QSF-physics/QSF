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
using Section = inipp::Ini<char>::Section;


#include "setup/config.h"
#include "absorber.h"
#include "coords.h"
#include "grid.h"
#include "wf.h"
#include "computations.h"
#include "dumps.h"
#include "average.h"
#include "hamiltonian.h"
#include "potential.h"
#include "field.h"
#include "coupling.h"
#include "propagator.h"
#include "routines.h"

// #include "setup/buffer.h"
namespace QSF
{
	void init(int argc, char* argv[])
	{
		MPI::init(argc, argv);
		logImportant("PROJECT: [%s] MPI PROCESSES: [%d]", STRINGIFY(PROJNAME), MPI::pSize);
		createDir(project_dir, results_dir);
		// createDir(home_dir, results_project_dir);
	}

	void finalize()
	{
		MPI_Finalize();
	}
};