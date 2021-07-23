#pragma once

// double* xcoords;			// x coords
// double* pcoords;			// p coords

// namespace Im
// {
// 	ind n;					// # of INITIAL STATE nodes (one direction),
// 	ind nn;					// [Derived] Square of Im::n
// 	ind m;				  	// [Derived] # of all ground state nodes (ex. Im::n*Im::n*Im::n for 3D)
// 	double L;				// Size of INITIAL STATE grid (one direction)
// };
// namespace Src
// {
// 	ind n;					// # of rows in prev MODE
// };

// double xmin;				  // [Derived] Lowest coordinate in grid = -L/2 (negative)
// double dx, dV, sqrt_dV;		  // [Derived] grid spacing, smallest volume and its square root
// double dt, dt2;				  // timestep (in atomic units) and half of timestep

// double A, B;				  // kinetic operator constants (gs)
// ind i, j, k;				  // shared temp indices
// //TEMPS
// ind readInd0, readInd1;		  // shared temp indices (MPI proceeses do not share memory)
// ind readInd2, readInd3;		  // shared temp indices (MPI proceeses do not share memory)
// double x, y, z;				  // temp coordinates vars
// ind ki, kj, kk, intemp;		  // temp local indices of op_ev_kin, op_ev_pot, energycalc

// double state_accuracy;		  // evolution stops when difference in energy becomes smaller than accuracy
// ind max_imaginary_steps;	  // maximal # of steps in IM
// int state;

// ind step;				  	  // for the main loop in both cases
// double timer = 0.;			  // Variable tracking real time in main.cpp
// ind ntsteps; 				  // steps in the main loop

// //OMP/MPI (local == current thread)


// ind m_l; 				// # of fftw_complex numbers to allocate for WF::psi
// ind n0_l;				// local # of rows
// ind n0_o;			// row # start
// ind n0_e;				// row # end, n0_e = n0_o + n0_l - 1;

// ind m_down, m_up;   		// locating min/max of initial state (psi0[])
// // double L;					// size of RE grid

// //RE GRID SIZE
// // ind n;						// full # of rows
// // ind m;						// full # of nodes
// // ind nn, n2;					// # rows squared, half # rows


// ind rowSize;				// 1 (1D), n (2D), nn (3D)
// double inv_m, sqrt_inv_m; 	// m variants
// double inv_dx;				// dx inverse
// double inv_2dx;				// (2*dx) inverse
// double dVm;					// dV/m

// //MOMENTA
// // double dp;					// Momentum space dx
// double dVP; 				// Momentum space volume element
// double sqrt_dVP, inv_dVP;	// dVP variants
// double sqrt_inv_dVP; 		// dVP variants: sqrt of inv dVP
// double p_scale;				// scales ψ(p,t) after FFTW on ψ(x,t)
// // double kin_scale;			// sets the kinetic energy scale
// double pmin;				// min momentum
// double pmax;				// max momentum
// cxd ctemp;				// shared temp kinetic and potential evolution operators

