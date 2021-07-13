#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)
// #include <sys/stat.h> //TODO: https://en.cppreference.com/w/cpp/filesystem/create_directory
#include <filesystem>

enum class IO_ATTR :char { READ = 'r', WRITE = 'w', APPEND = 'a' };
constexpr static char const* io_d[] = { "", "_x","_xy","_xyz" };
constexpr static char const* io_w[] = { "", "_bef","","","_aft" };
char parent_path[300];
char target_path[300];
char file_path[800];
char char_label[300];
auto home_dir = getenv("HOME");
auto scratch_dir = getenv("SCRATCH");
auto scratch_local_dir = getenv("SCRATCH_LOCAL");
auto project_dir = ".";
constexpr auto results_dir = "Results";
constexpr auto results_project_dir = "Results/" STRINGIFY(PROJNAME);
constexpr auto common_prefix = STRINGIFY(ELEC) "e" STRINGIFY(DIM) "d";
constexpr auto psi_ext = "psi";
constexpr auto log_ext = "log";
constexpr auto mom_ext = "mom";
constexpr auto pot_ext = "pot";
constexpr auto hhg_ext = "hhg";
constexpr auto flx_ext = "flx";
constexpr auto reg_ext = "reg";
constexpr auto exf_ext = "exf";
constexpr auto dat_ext = "dat";
constexpr auto eng_ext = "eng";

constexpr auto SEP = "_"; // Separator used in filenames


void createDir(char const* location, char const* subfolders)
{
	if (!MPI::rID)
	{
		sprintf(target_path, "%s/%s/", location, subfolders);
		// mkdir(target_path, 0777); //this won't work recursively
		std::filesystem::create_directories(target_path);
		// logIO("Created directory at %s", target_path);
	}
}
template <IO_ATTR ATTR>
FILE* fopen_with_check(char const* ffname)
{
	constexpr char atr[2] = { char(ATTR),0 };
	FILE* file_ptr = fopen(ffname, atr);
	if (file_ptr == nullptr)
	{
		logError("fopen_with_check: %s can not be opened (or more likely cannot be found!) with attributes %s", ffname, atr);
	}
	else { logIO("Opening [%s] file %s", atr, ffname); }
	return file_ptr;
}

template <DIMS D, typename WHEN>
void set_file_path(std::string_view name, char const* ext, bool binary)
{
	sprintf(file_path, "%s%s%s%s%s%s.%s%s%d",
			target_path,
			name.data(),
			common_prefix,
			SEP,
			WHEN::name,
			io_d[intDIMS(D)],
			ext,
			binary ? "b" : "",
			MPI::region);
}

template <typename WHEN, DIMS D, IO_ATTR ATTR>
FILE* openFile(std::string_view name, char const* ext, bool binary = false)
{
	if (!MPI::rID)
	{
		set_file_path<D, WHEN>(name, ext, binary);
		return fopen_with_check<ATTR>(file_path);
	}
	else return nullptr;
}
//Simplified
template <IO_ATTR ATTR>
FILE* openFile(std::string_view name, char const* ext, bool binary)
{
	return openFile <_TEMPORAL, DIMS(0), ATTR >(name, ext, binary);
}

template <IO_ATTR ATTR>
FILE* openOut(std::string_view name, const int index, bool binary)
{
	sprintf(char_label, "%s_%d", name.data(), index);
	return openFile<_TEMPORAL, DIMS(0), ATTR>(char_label, dat_ext, binary);
}

inline FILE* openLog(std::string_view name)
{
	return openFile<IO_ATTR::APPEND>(name, log_ext, false);
}

template <typename WHEN, REP R, DIMS D, IO_ATTR ATTR>
inline FILE* openPsi(std::string_view name, int state, ind step, bool binary)
{
	sprintf(char_label, "%s%d", name.data(), state);
	if (std::is_base_of_v<DURING, WHEN>) sprintf(char_label, "%s_%td", char_label, step);
	return openFile<WHEN, D, ATTR>(char_label,
								   R == REP::X ? psi_ext : mom_ext,
								   binary);
}

void closeFile(FILE* f) { if (!MPI::rID && f != nullptr) { fclose(f); } }

struct _Operator { OPTIMS maskOpt; };

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

template <typename T>
void dispatchComp(FILE* file, T t)
{
	char name[15];
	using types = decltype(t.types);
	if constexpr (0 < std::tuple_size_v<types>)
	{
		auto tup = t.types;
		ForEach(tup, [&](auto index)
				{
					char buf[15];
					auto sub = std::get<index>(tup);


					// sprintf(buf, t.name.data(), index.value);//, sub.name.data());
					// sprintf(name, "%14s", buf);
					// logInfo("Formatting %s", name);
					if constexpr (std::is_same_v<double, typename decltype(t)::type>)
					{
						if constexpr (std::is_base_of_v<_Operator, decltype(sub)>)
							sprintf(buf, "%s%s", t.name.data(), sub.name.data());
						if constexpr (!std::is_base_of_v<_Operator, decltype(sub)>)
							sprintf(buf, "%s%lu", t.name.data(), sub());
						sprintf(name, "%14s", buf);
						// logInfo("Formatting %s", name);
						fwrite(&name, sizeof(char), 15, file);
					}
					else if constexpr (std::is_same_v<cxd, typename decltype(t)::type>)
					{
						if constexpr (std::is_base_of_v<_Operator, decltype(sub)>)
							sprintf(buf, "%s%s%s", "RE_", t.name.data(), sub.name.data());
						if constexpr (!std::is_base_of_v<_Operator, decltype(sub)>)
							sprintf(buf, "%s%s%lu", "RE_", t.name.data(), sub());
						sprintf(name, "%14s", buf);
						// logInfo("Formatting %s", name);
						fwrite(&name, sizeof(char), 15, file);

						if constexpr (std::is_base_of_v<_Operator, decltype(sub)>)
							sprintf(buf, "%s%s%s", "IM_", t.name.data(), sub.name.data());
						if constexpr (!std::is_base_of_v<_Operator, decltype(sub)>)
							sprintf(buf, "%s%s%lu", "IM_", t.name.data(), sub());
						sprintf(name, "%14s", buf);
						// logInfo("Formatting %s", name);
						fwrite(&name, sizeof(char), 15, file);
					}
					else
					{
						strcpy(buf, name);
						constexpr size_t count = sizeof(typename decltype(t)::type) / sizeof(double);
						for (size_t i = 0; i < count; i++)
						{
							if constexpr (std::is_base_of_v<_Operator, decltype(sub)>)
								sprintf(buf, "#%zu%s%s", i, t.name.data(), sub.name.data());
							if constexpr (!std::is_base_of_v<_Operator, decltype(sub)>)
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


/** What the header contains
 * at the beginning: dim (# of dimensions written), n (# of points per dimension),
 * sizeof (size of the data, either sizeof(double) or sizeof(cxd),
 * min (real start of the axis, max (real end of the axis value),
 * delta (step between data in atomic units).
 * Total # of points expected is always pow(n,dim)
 */
template <DIMS D>
void writePsiBinaryHeader(FILE* file, double min, double max, double delta, DUMP_FORMAT F)
{
	if (file != nullptr)
	{
		// int dim = intDIMS(D);
		// fwrite(&dim, sizeof(int), 1, file);
		// fwrite(&n, sizeof(int), 1, file);
		// //If we start from dimension DIM, but compute distributions for
		// //dimension dim < DIM then the F should always be double (absolute square)
		// int size = (F.complex && (dim == DIM)) ? sizeof(cxd) : sizeof(double);
		// fwrite(&size, sizeof(int), 1, file);
		// fwrite(&min, sizeof(double), 1, file);
		// fwrite(&max, sizeof(double), 1, file);
		// fwrite(&delta, sizeof(double), 1, file);
	}
}

//This should match exactly the above writePsiBinaryHeader
bool readPsiBinaryHeader(FILE* file, ind expected_n)
{
	int dim;
	fread(&dim, sizeof(int), 1, file);
	// SKUBANY TEST POWODUJE EXC_BAD_ACCESS
	// logTest(dim == DIM, "Dimension of the input " psi_symbol" (%d) equals DIM (%d)", dim, DIM);
	int n_read;
	fread(&n_read, sizeof(int), 1, file);
	// logTest(n_read == expected_n, "1D Length of the input " psi_symbol " (%d) equals expected_n (%td)", n_read, expected_n);

	int size;
	fread(&size, sizeof(int), 1, file);
	//Discarded for now
	double tmp;
	fread(&tmp, sizeof(double), 1, file);
	fread(&tmp, sizeof(double), 1, file);
	fread(&tmp, sizeof(double), 1, file);
	return (sizeof(cxd) == size);
}

