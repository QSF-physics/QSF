#pragma once

#include <type_traits>

template<class T> constexpr inline T operator~(T a) { return (T) ~(uind)a; }
template<class T> constexpr inline T operator|(T a, T b) { return (T)((uind)a | (uind)b); }
template<class T> constexpr inline T operator&(T a, T b) { return (T)((uind)a & (uind)b); }
template<class T> constexpr inline T operator^(T a, T b) { return (T)((uind)a ^ (uind)b); }
template<class T> constexpr inline T& operator|=(T& a, T b) { return (T&)((uind&)a|= (uind)b); }
template<class T> constexpr inline T& operator&=(T& a, T b) { return (T&)((uind&)a&= (uind)b); }
template<class T> constexpr inline T& operator^=(T& a, T b) { return (T&)((uind&)a^= (uind)b); }
template<class T> constexpr inline T operator>>(T&& a, T&& b) { return (T)((uind)a >> (uind)b); }

/// @brief Mode of evolution during wf propagation
enum class MODE
{
	UNDEFINED= 0,
	IM			 = 1,
	RE			 = 2,
	ALL			 = IM | RE
};

#ifndef MODES_ENABLED
#define MODES_ENABLED MODE::ALL
#endif
/// @brief Compile-time filter which allows preparing single-task binaries
/// Currently used to prepare seperate imaginary- or real-time evolution binaries
#define MODE_FILTER_OPT(M) (bool(M & MODES_ENABLED))

constexpr DIMS operator""_D(unsigned long long val) { return static_cast<DIMS>(val); }

/// @brief Returns the name of a give propagation mode M
/// @param M (either IM or RE)
/// @return
static constexpr auto modeName(MODE M) { return M == MODE::IM ? "IM" : "RE"; }
enum class REP : DIMS
{
	NONE= 0,
	X		= 1,
	P		= 2,
	BOTH= 3
};

// constexpr REP inverse_rep(REP R) { return REP::BOTH ^ R; };

enum class AXIS : uind
{
	NO = 0,
	X	 = 1,
	Y	 = 1 << 1,
	Z	 = 1 << 2,	 // 3dim up to here
	// U = 1 << 3,
	// V = 1 << 4,
	// W = 1 << 5, //6dim up to here
	// R = 1 << 6,
	// S = 1 << 7,
	// T = 1 << 8, //9dim up to here
	XY = X | Y,
	XZ = X | Z,
	YZ = Y | Z,
	XYZ= X | Y | Z,
};

template<uind dir> constexpr AXIS Axis= AXIS(1 << dir);

constexpr AXIS getAxis(uind dir) { return AXIS(1 << dir); }

std::string axisName(AXIS ax)
{
	switch(ax)
	{
	case AXIS::NO: return "NONE";
	case AXIS::X: return "X";
	case AXIS::Y: return "Y";
	case AXIS::Z: return "Z";
	case AXIS::XY: return "XY";
	case AXIS::XZ: return "XZ";
	case AXIS::YZ: return "YZ";
	case AXIS::XYZ: return "XYZ";
	default: return std::to_string((uind)ax);
	}
}

enum class WHEN : ind
{
	AT_END	= -1,
	AT_START= 0,
	DURING	= 1,
};

// template <DIMS D>
// constexpr AXIS maxFreeAxis = AXIS((ind(1) << ind(D)) - 1);

// Returns number of directions used by AXIS
template<AXIS Ax>
static constexpr DIMS freeAxisCount=
	(uind(Ax) & 1) + (uind(Ax) == 0 ? 0 : freeAxisCount<(Ax >> AXIS(1))>);

struct _TEMPORAL
{
	static constexpr auto name= "";
};

template<typename... T> struct BEFORE: _TEMPORAL
{
	static constexpr auto name= "_bef";
};
template<typename... T> struct AFTER: _TEMPORAL
{
	static constexpr auto name= "_aft";
};

template<typename... T> struct DURING: _TEMPORAL
{};

template<size_t ST, typename... T> struct AT_STEP: _TEMPORAL
{
	static constexpr auto step= ST;
	// using types = tuple<T...>;
	static constexpr auto name= "_s";
};

template<size_t IN, typename... T> struct AT_INTERVAL: _TEMPORAL
{
	static constexpr auto interval= IN;
	// using types = tuple<T...>;
	static constexpr auto name		= "_i";
};

enum class OPTIMS
{
	NONE					= 0,
	PRECOMP_COORDS= BIT(1),

	PRECOMP_KIN					 = BIT(3),
	PRECOMP_VSTAT				 = BIT(4),
	PRECOMP_KIN_OP			 = BIT(5),
	PRECOMP_VSTAT_OP		 = BIT(6),
	VDYN_LINEAR_IN_X		 = BIT(7),
	VDYN_LINEAR_IN_Y		 = BIT(8),
	VDYN_LINEAR_IN_Z		 = BIT(9),
	VDYN_LINEAR_IN_ALL	 = VDYN_LINEAR_IN_X + VDYN_LINEAR_IN_Y + VDYN_LINEAR_IN_Z,
	STANDARD_OPTIMIATIONS= PRECOMP_KIN_OP + PRECOMP_VSTAT + PRECOMP_KIN +
												 PRECOMP_COORDS		// + VDYN_LINEAR_IN_ALL + PRECOMP_VSTAT_OP
		,
	VDYN_SUMXY_CONST	= BIT(10),
	VDYN_SUMXYZ_CONST = BIT(11),
	UNTRANSPOSED_P_DIM= BIT(12),
};

struct DUMP_FORMAT
{
	DIMS dim;
	REP rep				 = REP::X;
	bool complex	 = true;
	ind downscale	 = 1;
	bool evo_backup= false;
};
// //Return the flag for current DIM
// // constexpr DIMS DIMflag = (DIM == 1 ? D1 : (DIM == 2 ? D2 : D3));
// constexpr int intDIMS(DIMS D)
// {
// 	int v = int(D);
// 	int MultiplyDeBruijnBitPosition[32] =
// 	{
// 	  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
// 	  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
// 	};
// 	return 1 + MultiplyDeBruijnBitPosition[((uint32_t)((v & -v) * 0x077CB531U)) >> 27];
// 	   // return D == DIMS::D1 ? 1 : (D == DIMS::D2 ? 2 : D == DIMS::D3 ? 3 : 0);
// }
// enum class DIMS
// {
// 	D0 = 0,
// 	D1 = 1,
// 	D2 = 2,
// 	D3 = 3,
// 	D4 = 4,
// 	D5 = 5,
// 	D6 = 6,
// 	D7 = 7,
// 	D8 = 8,
// 	D9 = 9,
// 	ALL = D1 + D2 + D3
// };
/* This allows for filtering  */
// #ifndef ONLY_DIM
// #define ONLY_DIM DIMS::ALL
// #endif

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