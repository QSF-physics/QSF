#include <fstream>

namespace QSF
{
	using Section = inipp::Ini<char>::Section;

	struct Config
	{
		std::string source_section;
		inipp::Ini<char> ini;
		char config[100];

		Config() = default;
		explicit Config(std::string section, std::string_view name)
		{
			logInfo("Parsing %s", name.data());
			source_section = section;
			// if (MPI::pID)
			std::ifstream is(IO::project_dir / std::string(name));
			// ini.clear();
			ini.parse(is);
			// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
			ini.strip_trailing_comments();
			ini.default_section(ini.sections["DEFAULT"]);
			ini.interpolate();
			logINI("Parsed & interpolated project.ini file:");
			ini.generate(std::cout);

			// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
			// if (DEBUG & DEBUG_INI) MPI_Barrier(MPI_COMM_WORLD);
		}
		void subdirectory(std::filesystem::path sub)
		{
			if (!MPI::rID)
				std::filesystem::create_directories(std::filesystem::current_path() / sub);
			std::filesystem::current_path(std::filesystem::current_path() / sub);
		}
	};
};