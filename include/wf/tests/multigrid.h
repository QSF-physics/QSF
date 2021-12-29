#pragma once
#ifdef TEST_MULTIGRID_MASK
FILE* mask_debug = nullptr;
bool mask_debug_done = false;

FILE* slice_mask_debug = nullptr;
bool slice_mask_debug_done = false;


#define TEST_MULTIGRID_MASK_P0 if (mask_debug != nullptr) { fflush(mask_debug); fclose(mask_debug); mask_debug = nullptr; mask_debug_done = true; }\
	if (!mask_debug_done) mask_debug = IO::fopen_with_check("multigrid_mask_region" + std::to_string(MPI::region) + "_rID" + std::to_string(MPI::rID) + ".csv", "w"); \
	if (MPI::group == 1) {\
	if (slice_mask_debug != nullptr) { fflush(slice_mask_debug); fclose(slice_mask_debug); slice_mask_debug = nullptr; slice_mask_debug_done = true; }\
	if (!slice_mask_debug_done) slice_mask_debug = IO::fopen_with_check("multigrid_slice_mask_region" + std::to_string(MPI::region) + "_rID" + std::to_string(MPI::rID) + ".csv", "w");}

#define TEST_MULTIGRID_MASK_P1 if (!mask_debug_done) { \
	fprintf(mask_debug, DIM == 3 ? "%td, %td, %td, %g\n" : "%td, %td, %g\n", (counters[dirs] + (dirs ? 0 : pos_lx.first))..., region_mask); }

#define TEST_MULTIGRID_SLICE_MASK_P1 if (MPI::group == 1){ if (!slice_mask_debug_done) { \
	fprintf(slice_mask_debug, DIM == 3 ? "%td, %td, %td, %g, %g\n" : "%td, %td, %g, %g\n", (counters[dirs] + (dirs ? 0 : pos_lx.first) + slider<dirFree, dirs>(boxIndex))..., mask_value, mask_value_corr); }}

#define TEST_MULTIGRID_SLICE_MASK_P2 if (MPI::group == 1) { if (!slice_mask_debug_done) { \
fprintf(slice_mask_debug, DIM == 3 ? "%td, %td, %td, %g, %g\n" : "%td, %td, %g, %g\n", (counters[dirs] + (dirs ? 0 : pos_lx.first) + tslider<dirFree, dirs>(boxIndex))..., mask_value, mask_value_corr); }}

#else
#define TEST_MULTIGRID_MASK_P0 
#define TEST_MULTIGRID_MASK_P1 
#define TEST_MULTIGRID_SLICE_MASK_P1
#define TEST_MULTIGRID_SLICE_MASK_P2
#endif