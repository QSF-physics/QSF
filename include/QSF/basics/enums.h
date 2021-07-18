#pragma once

#include <type_traits>

template<class T> constexpr inline T operator~ (T a) { return (T)~(unsigned)a; }
template<class T> constexpr inline T operator| (T a, T b) { return (T)((unsigned)a | (unsigned)b); }
template<class T> constexpr inline T operator& (T a, T b) { return (T)((unsigned)a & (unsigned)b); }
template<class T> constexpr inline T operator^ (T a, T b) { return (T)((unsigned)a ^ (unsigned)b); }
template<class T> constexpr inline T& operator|= (T& a, T b) { return (T&)((unsigned&)a |= (unsigned)b); }
template<class T> constexpr inline T& operator&= (T& a, T b) { return (T&)((unsigned&)a &= (unsigned)b); }
template<class T> constexpr inline T& operator^= (T& a, T b) { return (T&)((unsigned&)a ^= (unsigned)b); }

struct DUMP_FORMAT
{
	bool binary = true;
	bool complex = true;
	bool unnormalized = false;
	bool no_coords = true;
	bool subtract_initial_wf = false;
};

struct DATA_FORMAT
{
	bool binary = true;
	// constexpr DATA_FORMAT(bool binary = false) : binary(binary) {}
};

enum class MODE
{
	UNDEFINED = 0,
	IM = 1,
	RE = 2,
	ALL = IM + RE
};

#ifndef ONLY_MODE
#define ONLY_MODE MODE::ALL
#endif

enum class DIMS
{
	D0 = 0,
	D1 = 1,
	D2 = 2,
	D3 = 4,
	ALL = D1 + D2 + D3
};

/* This allows for filtering  */
#ifndef ONLY_DIM
#define ONLY_DIM DIMS::ALL
#endif

// enum class ELEC
// {
// }


// template <MODE M> constexpr string_view modeName()
// {
// 	if constexpr (M == IM) return "IM";
// 	if constexpr (M == RE) return "RE";
// 	if constexpr (M == EX) return "EX";
// 	return "";
// }
static constexpr auto modeName(MODE M) { return M == MODE::IM ? "IM" : "RE"; }
enum class REP
{
	NONE = 0,
	X = 1,
	P = 2,
	BOTH = 3
};

constexpr REP invREP(REP R) { return REP::BOTH ^ R; };


enum class AXIS
{
	NO = 0,
	X = 1,
	Y = 2,
	Z = 4,
	XY = X | Y,
	XZ = X | Z,
	YZ = Y | Z,
	XYZ = X | Y | Z,
// #if DIM ==1
// 	ALL = X,
// 	LAST = X
// #elif DIM == 2
// 	ALL = XY,
// 	LAST = Y
// #elif DIM == 3
// 	ALL = XYZ,
// 	LAST = Z
// #endif
};

enum class WHEN : ind
{
	AT_START = 0,
	AT_END = -1,
	DURING,
};

using FREE_COORD = AXIS;

template <DIMS D>
constexpr FREE_COORD maxFreeCoord = FREE_COORD((1 << int(D)) - 1);

// Returns number of dimensions used by direction (size_t)
template <FREE_COORD fc>
static constexpr int freeCoordCount = (fc == FREE_COORD::NO ? 0 :
	((fc == FREE_COORD::X | fc == FREE_COORD::Y | fc == FREE_COORD::Z) ? 1 :
	 (fc == FREE_COORD::XYZ ? 3 : 2)));


// template <FREE_COORD...FC>
// struct FreeCoords
// {
// 	static constexpr auto size = sizeof...(FC);
// }
// static int freeCoordDir(FREE_COORD fc)
// {
// 	return (int(fc) == 4) ? 2 : ((int)(fc)-1);
// }
// Returns number of dimensions used by direction (DIMS flag)



struct _TEMPORAL { static constexpr auto name = ""; };

template <typename...T>
struct BEFORE : _TEMPORAL
{
	static constexpr auto name = "_bef";
};
template <typename...T>
struct AFTER : _TEMPORAL
{
	static constexpr auto name = "_aft";
};

template <typename ... T>
struct DURING : _TEMPORAL
{

};

template <size_t ST, typename ... T>
struct AT_STEP : _TEMPORAL
{
	static constexpr auto step = ST;
	// using types = tuple<T...>;
	static constexpr auto name = "_s";
};

template <size_t IN, typename ... T>
struct AT_INTERVAL : _TEMPORAL
{
	static constexpr auto interval = IN;
	// using types = tuple<T...>;
	static constexpr auto name = "_i";
};

struct EARLY
{

	template <typename WHEN> constexpr bool canRun()  const noexcept { return std::is_same_v<WHEN, EARLY>; }
};
struct LATE
{
	template <typename WHEN> constexpr bool canRun()  const noexcept { return std::is_same_v<WHEN, LATE>; }
};


//Return the flag for current DIM
// constexpr DIMS DIMflag = (DIM == 1 ? D1 : (DIM == 2 ? D2 : D3));
constexpr int intDIMS(DIMS D) { return D == DIMS::D1 ? 1 : (D == DIMS::D2 ? 2 : D == DIMS::D3 ? 3 : 0); }


enum class OPTIMS
{
	NONE = 0,
	PRECOMP_COORDS = BIT(1),

	PRECOMP_KIN = BIT(3),
	PRECOMP_VSTAT = BIT(4),
	PRECOMP_KIN_OP = BIT(5),
	PRECOMP_VSTAT_OP = BIT(6),
	VDYN_LINEAR_IN_X = BIT(7),
	VDYN_LINEAR_IN_Y = BIT(8),
	VDYN_LINEAR_IN_Z = BIT(9),
	VDYN_LINEAR_IN_ALL = VDYN_LINEAR_IN_X + VDYN_LINEAR_IN_Y + VDYN_LINEAR_IN_Z,
	STANDARD_OPTIMIATIONS = PRECOMP_KIN_OP + PRECOMP_VSTAT + PRECOMP_KIN + PRECOMP_COORDS// + VDYN_LINEAR_IN_ALL + PRECOMP_VSTAT_OP
	, VDYN_SUMXY_CONST = BIT(10),
	VDYN_SUMXYZ_CONST = BIT(11),
	UNTRANSPOSED_P_DIM = BIT(12),
};
