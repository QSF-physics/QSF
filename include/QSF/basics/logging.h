#define DEBUG_CALLS BIT(0) 
#define DEBUG_TESTS BIT(1)
#define DEBUG_SETUP BIT(2)
#define DEBUG_ALLOC BIT(3)
#define DEBUG_IO BIT(4)
#define DEBUG_INI BIT(5)
#define DEBUG_MPI BIT(6)
#define DEBUG_BUFFER BIT(7)
#define DEBUG_DUMPS BIT(8)

#define GOBACK "" //"../../"
#define FAIL() {  MPI_Barrier(MPI_COMM_WORLD); MPI_Finalize();  exit(EXIT_FAILURE); }
FILE* file_log;

// void* operator new (size_t size, char const* filename, int line) {
// 	void* ptr = new char[size];
// 	cout << "size = " << size << " filename = " << filename << " line = " << line << endl;
// 	return ptr;
// }

// #define new new(__FILE__, __LINE__)

#if defined DEBUG && DEBUG
int qdd;
#define qd() {MPI_Barrier(MPI_COMM_WORLD); if (!MPI::pID) { fprintf(stdout, "%27s:%-5i| ", GOBACK __FILE__, __LINE__);  fprintf(stdout, "------------------------------------------------------------------------------ >QD: %d\n", qdd); qdd++; }}
// use __func__

#else
#define qd() {}
#endif

#if defined _WIN32 || defined _LOG_NO_COLOR
#define __LOG_RED ""
#define __LOG_GREEN ""
#define __LOG_BLUE ""
#define __LOG_MAGENTA ""
#define __LOG_CYAN ""
#define __LOG_YELLOW ""
#define __LOG_WHITE ""
#define __LOG_NC ""
#else
#define __LOG_RED "\033[1;31m"
#define __LOG_GREEN "\033[0;32m"
#define __LOG_BLUE "\033[1;34m"
#define __LOG_MAGENTA "\033[1;35m"
#define __LOG_CYAN "\033[1;36m"
#define __LOG_YELLOW "\033[0;33m"
#define __LOG_WHITE "\033[0;37m"
#define __LOG_NC "\033[m"
#endif
#define __LOG_NL "\n"

//temp hack
constexpr const char* file_name(const char* path)
{
	const char* file = path + 43;
	// while (*path) {
	// 	if (*path++ == '/') {
	// 		file = path;
	// 	}
	// }
	return file;
}


// WORD ABOUT THE NOTION: __ performs no checks, _ performs file check, but no MPI::pID check
#ifdef _NSL
#define __LOG2S(format, ...) {}
#else
#define __LOG2S(format, ...) { printf(format, ##__VA_ARGS__); }
#endif
#ifdef _NFL
#define __LOG2F(...) {}
#else
#define __LOG2F(...) {fprintf(file_log, __VA_ARGS__);}
#endif

#if (DEBUG & DEBUG_CALLS)
#define __LOG_SRC_LOC(color) {__LOG2S("%31s:%-4i| %s" , file_name(__FILE__), __LINE__, color); }
#else
#define __LOG_SRC_LOC(color) {__LOG2S("%s", color); }
// #define __LOG_SRC_LOC(color) {}
#endif 

#define _LOG_INLINE_START(color) { __LOG_SRC_LOC(color) }
#define LOG_INLINE_START(color) { if (!MPI::pID) { _LOG_INLINE_START(color) } }
#define __LOG_INLINE_END() { __LOG2S(__LOG_NC __LOG_NL) __LOG2F(__LOG_NL) }
#define _LOG_INLINE_END() { __LOG2S(__LOG_NC __LOG_NL)  if (file_log != nullptr) __LOG2F(__LOG_NL) }
#define LOG_INLINE_END() {if (!MPI::pID) { _LOG_INLINE_END() }}

#define __LOG_INLINE(format, ...) { __LOG2S(format, ##__VA_ARGS__) __LOG2F(format, ##__VA_ARGS__) }
#define _LOG_INLINE(format, ...) { __LOG2S(format, ##__VA_ARGS__) if (file_log != nullptr) __LOG2F(format, ##__VA_ARGS__) }
#define LOG_INLINE(format, ...) { if (!MPI::pID) { _LOG_INLINE(format, ##__VA_ARGS__) }}

#define _logColor(color, format, ...) {_LOG_INLINE_START(color) _LOG_INLINE(format, ##__VA_ARGS__) _LOG_INLINE_END()}
#define logColor(color, format, ...) {if (!MPI::pID) { _logColor(color, format, ##__VA_ARGS__) }}

#define logInfo(format, ...) { logColor("", format, ##__VA_ARGS__)}
#define logUser(format, ...) { logColor(__LOG_CYAN, format, ##__VA_ARGS__)}
#define logImportant(format, ...) { logColor(__LOG_BLUE, format, ##__VA_ARGS__)}
#define logWarning(format, ...) { logColor(__LOG_MAGENTA, format, ##__VA_ARGS__) }

#define __TS(res) (bool(res) ? check_symbol : "x")
#define __TC(res) (bool(res) ? __LOG_GREEN : __LOG_RED)
#define __TC2(res) (bool(res) ? __LOG_WHITE : __LOG_RED)

#define plogTest(res, ...) {if (!MPI::pID) { _LOG_INLINE_START(__TC2(res)) _LOG_INLINE(__VA_ARGS__) __LOG2F(": [%s]", __TS(res)) __LOG2S(": [%s%s%s]", __TC(res), __TS(res), __TC2(res)) _LOG_INLINE_END() }}
#define plogTestFatal(res, ...) { plogTest(res,__VA_ARGS__) if(! bool(res)) { FAIL() } }

#define logTest(res, ...) {if ((DEBUG & DEBUG_TESTS)|| ! bool(res)) {plogTest(res, __VA_ARGS__) }}
#define logTestFatal(res, ...) {if ((DEBUG & DEBUG_TESTS)|| ! bool(res)) {plogTestFatal(res, __VA_ARGS__) }}

#define logSuccess(format, ...) { logColor(__LOG_GREEN, format, ##__VA_ARGS__)}
#define logError(...) { logColor(__LOG_RED, __VA_ARGS__) FAIL() }



#define logSETUP(...) {if constexpr (bool(DEBUG & DEBUG_SETUP)) { logColor(__LOG_BLUE, __VA_ARGS__) }}
#define logALLOC(...) {if constexpr (bool(DEBUG & DEBUG_ALLOC)) { logColor(__LOG_RED, __VA_ARGS__) }}
#define logIO(...) {if constexpr (bool(DEBUG & DEBUG_IO)) { logColor(__LOG_YELLOW, __VA_ARGS__) }}
#define logINI(...) {if constexpr (bool(DEBUG & DEBUG_INI)) { logColor(__LOG_WHITE, __VA_ARGS__) }}

#define logMPI(...) { if constexpr (bool(DEBUG & DEBUG_MPI)) {  logColor(__LOG_YELLOW, __VA_ARGS__)   }}
#define logDUMPS(...) {if constexpr (bool(DEBUG & DEBUG_DUMPS)) { logColor(__LOG_CYAN, __VA_ARGS__) }}
#define logBUFFER(...) {if constexpr (bool(DEBUG & DEBUG_BUFFER)) { logColor(__LOG_YELLOW, __VA_ARGS__) }}
