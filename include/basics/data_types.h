#pragma once
#include <complex>

/* I hope that what follows will not add much confusion:
here are only 3 frequently using basic data types */

/* Used throughout the project as array datatype  */
using cxd = std::complex<double>;
/* Array indexing proper for large array and FFTW */
using ind = ptrdiff_t;
/* Unsigned version of the above  */
using uind = size_t;
/* "256 dimensions should be enough for everybody" - Bill Gates, 2021 */
using DIMS = unsigned char;
