
#define STRINGIFY_(x) #x
#define STRINGIFY(x)	STRINGIFY_(x)

#include <cstring>
struct _Operator
{
	OPTIMS maskOpt;
};

namespace QSF
{
	namespace IO
	{

		using namespace std::filesystem;
		using namespace std::string_literals;
		path root_dir;
		const path home_dir{std::getenv("HOME")};
		const path project_dir= current_path();
		const path project_name{STRINGIFY(PROJNAME)};
		const path results_dir{"Results"};
		const path temp_dir						= temp_directory_path();
		const path results_project_dir= results_dir / project_name;

		constexpr auto SEP		= "_";
		constexpr auto psi_ext= ".psi";
		constexpr auto log_ext= ".log";
		constexpr auto mom_ext= ".mom";
		constexpr auto pot_ext= ".pot";
		constexpr auto hhg_ext= ".hhg";
		constexpr auto flx_ext= ".flx";
		constexpr auto reg_ext= ".reg";
		constexpr auto exf_ext= ".exf";
		constexpr auto dat_ext= ".dat";
		constexpr auto eng_ext= ".eng";

		FILE* fopen_with_check(std::string ffname, const char* atr)
		{

			FILE* file_ptr= fopen(ffname.c_str(), atr);
			if(file_ptr == nullptr)
			{
				logError("fopen_with_check: %s can not be opened (or more likely cannot be found!) with "
								 "attributes %s",
								 ffname.c_str(),
								 atr);
			}
			else { logIO("Opening [%s] file %s", atr, ffname.c_str()); }
			return file_ptr;
		}

		// FILE* openOut(std::string_view name, bool binary)
		// {
		// 	sprintf(char_label, "%s_%d", name.data(), index);
		// 	return fopen(ffname, atr);

		// 	openFile<_TEMPORAL, DIMS(0), ATTR>(char_label, dat_ext, binary);
		// }

		// enum class IO_ATTR :char { READ = 'r', WRITE = 'w', APPEND = 'a' };
		// constexpr static char const* io_d[] = { "", "_x","_xy","_xyz" };
		// constexpr static char const* io_w[] = { "", "_bef","","","_aft" };
		// char parent_path[300];

		// char file_path[800];
		// char char_label[300];
		// constexpr auto common_prefix = "";//STRINGIFY(ELEC) "e" STRINGIFY(DIM) "d";

		// template <DIMS D, typename WHEN>
		// void set_file_path(std::string_view name, char const* ext, bool binary)
		// {
		// 	sprintf(file_path, "%s/%s%s%s%s%s.%s%s%d",
		// 			target_path.c_str(),
		// 			name.data(),
		// 			common_prefix,
		// 			SEP,
		// 			WHEN::name,
		// 			io_d[D],
		// 			ext,
		// 			binary ? "b" : "",
		// 			MPI::region);
		// }

		// template <typename WHEN, DIMS D, IO_ATTR ATTR>
		// FILE* openFile(std::string_view name, char const* ext, bool binary = false)
		// {
		// 	if (!MPI::rID)
		// 	{
		// 		set_file_path<D, WHEN>(name, ext, binary);
		// 		return fopen_with_check<ATTR>(file_path);
		// 	}
		// 	else return nullptr;
		// }
		// //Simplified
		// template <IO_ATTR ATTR>
		// FILE* openFile(std::string_view name, char const* ext, bool binary)
		// {
		// 	return openFile <_TEMPORAL, DIMS(0), ATTR >(name, ext, binary);
		// }

		// template <IO_ATTR ATTR>
		// FILE* openOut(std::string_view name, const int index, bool binary)
		// {
		// 	sprintf(char_label, "%s_%d", name.data(), index);
		// 	return openFile<_TEMPORAL, DIMS(0), ATTR>(char_label, dat_ext, binary);
		// }

		// inline FILE* openLog(std::string_view name)
		// {
		// 	return openFile<IO_ATTR::APPEND>(name, log_ext, false);
		// }

		// template <typename WHEN, REP R, DIMS D, IO_ATTR ATTR>
		// inline FILE* openPsi(std::string_view name, int state, ind step, bool binary)
		// {
		// 	sprintf(char_label, "%s%d", name.data(), state);
		// 	if (std::is_base_of_v<DURING<>, WHEN>) sprintf(char_label, "%s_%td", char_label, step);
		// 	return openFile<WHEN, D, ATTR>(char_label,
		// 								   R == REP::X ? psi_ext : mom_ext,
		// 								   binary);
		// }

		// void closeFile(FILE* f) { if (!MPI::rID && f != nullptr) { fclose(f); } }

		// template <typename A, typename B>
		// inline constexpr std::string_view buildName(A a, B b)
		// {

		// 	if constexpr (std::is_base_of_v<_Operator, B>) return join_v<a.name, b.name>;
		// 	if constexpr (!std::is_base_of_v<_Operator, B>)
		// 	{
		// 		auto c = std::string_view{ b() };
		// 		return join_v < a.name, to_string(c)>;
		// 	}
		// }

		template<typename T> void dispatchComp(FILE* file, T t)
		{
			char name[15];
			using types= decltype(t.types);
			if constexpr(0 < std::tuple_size_v<types>)
			{
				auto tup= t.types;
				ForEach(tup,
								[&](auto index)
								{
									char buf[15];
									auto sub= std::get<index>(tup);

									// sprintf(buf, t.name.data(), index.value);//, sub.name.data());
									// sprintf(name, "%14s", buf);
									// logInfo("Formatting %s", name);
									if constexpr(std::is_same_v<double, typename decltype(t)::type>)
									{
										if constexpr(std::is_base_of_v<_Operator, decltype(sub)>)
											sprintf(buf, "%s%s", t.name.data(), sub.name.data());
										if constexpr(!std::is_base_of_v<_Operator, decltype(sub)>)
											sprintf(buf, "%s%lu", t.name.data(), sub());
										sprintf(name, "%14s", buf);
										// logInfo("Formatting %s", name);
										fwrite(&name, sizeof(char), 15, file);
									}
									else if constexpr(std::is_same_v<cxd, typename decltype(t)::type>)
									{
										if constexpr(std::is_base_of_v<_Operator, decltype(sub)>)
											sprintf(buf, "%s%s%s", "RE_", t.name.data(), sub.name.data());
										if constexpr(!std::is_base_of_v<_Operator, decltype(sub)>)
											sprintf(buf, "%s%s%lu", "RE_", t.name.data(), sub());
										sprintf(name, "%14s", buf);
										// logInfo("Formatting %s", name);
										fwrite(&name, sizeof(char), 15, file);

										if constexpr(std::is_base_of_v<_Operator, decltype(sub)>)
											sprintf(buf, "%s%s%s", "IM_", t.name.data(), sub.name.data());
										if constexpr(!std::is_base_of_v<_Operator, decltype(sub)>)
											sprintf(buf, "%s%s%lu", "IM_", t.name.data(), sub());
										sprintf(name, "%14s", buf);
										// logInfo("Formatting %s", name);
										fwrite(&name, sizeof(char), 15, file);
									}
									else
									{
										strcpy(buf, name);
										constexpr size_t count= sizeof(typename decltype(t)::type) / sizeof(double);
										for(size_t i= 0; i < count; i++)
										{
											if constexpr(std::is_base_of_v<_Operator, decltype(sub)>)
												sprintf(buf, "#%zu%s%s", i, t.name.data(), sub.name.data());
											if constexpr(!std::is_base_of_v<_Operator, decltype(sub)>)
												sprintf(buf, "#%zu%s%lu", i, t.name.data(), sub());
											sprintf(name, "%14s", buf);
											// logInfo("Formatting %s", name);
											fwrite(&name, sizeof(char), 15, file);
										}
									}
								});
			}
			// else
			// {
			// 	// sprintf(name, "%14s", t.name.data());
			// 	// logInfo("Formatting single %s", name);
			// 	// // logInfo("Formatting %s", name);
			// 	// fwrite(&name, sizeof(char), 15, file);
			// }
		}

		void writePsiBinaryHeader(FILE* file, ind* ns, double* dxs, bool* bounded, DUMP_FORMAT df)
		{
			if(file != nullptr)
			{
				// fwrite(&df.dim, sizeof(DIMS), 1, file);
				// fwrite(&df.rep, sizeof(DIMS), 1, file);
				// fwrite(&df.unnormalized, sizeof(bool), 1, file);
				// fwrite(&df.initial_wf_subtracted, sizeof(bool), 1, file);
				// //If we start from dimension DIM, but compute distributions for
				// //dimension dim < DIM then the F should always be double (absolute square)
				// int size = true ? sizeof(cxd) : sizeof(double);
				// fwrite(&size, sizeof(int), 1, file);

				// fwrite(ns, sizeof(ind), df.dim, file);
				// fwrite(dxs, sizeof(double), df.dim, file);
				// fwrite(bounded, sizeof(bool), df.dim, file);
			}
		}

		// This should match exactly the above writePsiBinaryHeader
		template<DIMS DIM> bool readPsiBinaryHeader(FILE* file)
		{
			DIMS dim, rep;
			fread(&dim, sizeof(DIMS), 1, file);
			fread(&rep, sizeof(DIMS), 1, file);
			bool unnormalized, initial_wf_subtracted;
			fread(&unnormalized, sizeof(bool), 1, file);
			fread(&initial_wf_subtracted, sizeof(bool), 1, file);
			int size;
			fread(&size, sizeof(int), 1, file);
			ind ns[DIM];
			double dxs[DIM];
			bool bounded[DIM];
			fread(&ns, sizeof(ind), DIM, file);
			fread(&dxs, sizeof(double), DIM, file);
			fread(&bounded, sizeof(bool), DIM, file);
			// SKUBANY TEST POWODUJE EXC_BAD_ACCESS
			// logTest(dim == DIM, "Dimension of the input " psi_symbol" (%d) equals DIM (%d)", dim, DIM);

			// logTest(n_read == expected_n, "1D Length of the input " psi_symbol " (%d) equals expected_n
			// (%td)", n_read, expected_n);

			// int size;
			// fread(&size, sizeof(int), 1, file);
			// Discarded for now
			// double tmp;
			// fread(&tmp, sizeof(double), 1, file);
			// fread(&tmp, sizeof(double), 1, file);
			// fread(&tmp, sizeof(double), 1, file);
			return (sizeof(cxd) == size);
		}
	};	 // namespace IO
};		 // namespace QSF