#pragma once

#ifdef NO_INPUT_FILE
#define constex constexpr
#else
#define NO_INPUT_FILE 0
#define constex
#endif 


// FORMAT CONSTANTS
#define STATE_NAME_PADDING ".2"
#define FMT_TIME "%15.4g"
#define FMTS_TIME "%15s"
#define FMT_D "%.15g"
#define FMT_DOUBLE "|%13.13g "
#define FMTS_DOUBLE "%15s"
#define FMT_STEPS "%7td/%7td"
#define FMT_STEPS_IM "%15td"
#define FMT_STEP "%7td"
#define FMT_ETA "%9d:%02d:%02d"
#define FMTS_ETA "%15s"
#define FMTS_AMPLITUDES "%13s%2d|"
#define FMTS_AVG "%13s %1s"

#define FMT_SEP "\t"
#define FMT_END "\n"
#define FMTS_POS_NAME "%8s | "
#define FMTS_POS_NAME_ITER "%6s %-1s | "
#define FMTS_POS_NAME_ITER_AMPL " %5s%2d| | "
#define FMTS_POS_NUM "%8d | "
#define psi_symbol "ψ"
#define check_symbol "✓"

//TODO: modernize
#define unit_intensity  3.50944e16 //Unit intensity
constexpr cxd I = { 0.0, 1.0 };	// Imaginary unit
constexpr auto eV = 27.2113; // ElectonVolt unit conversion factor
constexpr double pi = M_PI;
constexpr auto femtos = 41.34138;	  		// Femtosecond is this many atomic units


constexpr double pi2 = pi / 2.0;
constexpr double inv_twopi = 1.0 / pi / 2.0;
constexpr double inv_pi = 1.0 / pi;
constexpr double twopi = 2.0 * pi;
constexpr double threepi2 = 1.5 * pi;