inipp::Ini<char> ini;
#include <fstream>

void useProjectIni(std::string_view name)
{
	// logInfo("Parsing %s", name.data());
	// if (MPI::pID)
	std::ifstream is(name.data());
	ini.parse(is);
	// logINI("Raw INI file:");
	// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
	ini.strip_trailing_comments();
	constexpr auto dimelec = STRINGIFY(ELEC) "e" STRINGIFY(DIM) "d";
	ini.default_section(ini.sections[dimelec]);
	ini.default_section(ini.sections["DEFAULT"]);
	ini.interpolate();
	// logINI("Parsed & interpolated project.ini file:");
	// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
	// if (DEBUG & DEBUG_INI) MPI_Barrier(MPI_COMM_WORLD);
}

