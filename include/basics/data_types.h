#pragma once
#include <complex>

/* I hope that what follows will not add much confusion:
here are only 3 frequently using basic data types */

/// @brief Complex double shorthand used throughout the project as array datatype
using cxd = std::complex<double>;
/// @brief Array indexing proper for large array and FFTW
using ind = ptrdiff_t;
/// @brief Unsigned version of ind (large array index)
using uind= size_t;

/// @brief Small dimensional index
/// - "256 dimensions should be enough for everybody" - Bill Gates, 2021
using DIMS= unsigned char;
