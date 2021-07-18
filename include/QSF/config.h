#include <fstream>
struct Config
{
	inipp::Ini<char> ini;
	char config[100];

	Config() = default;
	explicit Config(std::string_view name, int DIMS, int ELEC)
	{
		logInfo("Parsing %s", name.data());
		// if (MPI::pID)
		std::ifstream is(name.data());
		// ini.clear();
		ini.parse(is);
		// logINI("Raw INI file:");
		// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
		ini.strip_trailing_comments();
		sprintf(config, "%de%dd", DIMS, ELEC);
		// constexpr auto dimelec = STRINGIFY(ELEC) "e" STRINGIFY(DIM) "d";
		ini.default_section(ini.sections[config]);
		ini.default_section(ini.sections["DEFAULT"]);
		ini.interpolate();
		// logINI("Parsed & interpolated project.ini file:");
		// ini.generate(std::cout);
		// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
		// if (DEBUG & DEBUG_INI) MPI_Barrier(MPI_COMM_WORLD);
	}
};