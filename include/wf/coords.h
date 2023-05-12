#pragma once
namespace QSF
{
	template<DIMS dims> struct CoordinateSystem
	{
		double dx[dims];
		ind n[dims];
		static constexpr DIMS DIM= dims;
	};

	/// @brief Generic N-dimensional Cartesian grid
	/// @tparam MPIRegions
	/// @tparam dims
	template<DIMS dims, class MPIRegions> struct CartesianGrid_: CoordinateSystem<dims>
	{
		// static_assert(std::is_base_of_v< AbsorberType, AT>,
		// "Second CartesianGrid template argument must be type derived from AbsorberType");

		using MPIGridComm= MPIRegions;
		using MPIDivision= typename MPIRegions::MPIDivision;
		/* Init of the above part is postponed to the wf */

		using CoordinateSystem<dims>::n;
		using CoordinateSystem<dims>::dx;
		using CoordinateSystem<dims>::DIM;
		// Read n/2, holds counts of half-grid nodes
		ind n2[dims];
		// Read n*n, holds squares counts (n) of grid nodes
		ind nn[dims];
		// Total count of grid nodes
		ind m;

		double inv_m;
		double dV;
		double dVm;
		double L[dims];
		double xmin[dims];
		double inv_dx[dims];
		double inv_2dx[dims];
		double inv_n[dims];

		double dp[dims];
		double pmin[dims];
		double pmax[dims];
		double kin_scale[dims];
		double dVP;

		template<REP R> double vol()
		{
			if constexpr(R == REP::X) return dV;
			else
				return dVm;
		}

		template<uind... Is> void init(seq<Is...>)
		{
			m= 1;
			(
				[&]
				{
					n2[Is]	 = n[Is] / 2;
					nn[Is]	 = n[Is] * n[Is];
					inv_n[Is]= 1.0 / n[Is];
					m*= n[Is];
				}(),
				...);
			inv_m= 1.0 / m;

			dV = 1.0;
			dVP= 1.0;
			(
				[&]
				{
					dV*= dx[Is];
					L[Is]				 = (dx[Is] * (n[Is] - 1));
					xmin[Is]		 = -0.5 * L[Is];
					inv_dx[Is]	 = 1.0 / dx[Is];
					inv_2dx[Is]	 = inv_dx[Is] / 2.0;
					dp[Is]			 = 2.0 * pi / double(n[Is]) / dx[Is];
					pmin[Is]		 = -pi / dx[Is];
					pmax[Is]		 = -pmin[Is] - dp[Is];
					kin_scale[Is]= dp[Is] * dp[Is] * 0.5;
					dVP*= dp[Is];
				}(),
				...);
			dVm= dV * inv_m;

			// logSETUP("CartesianGrid init");
			logSETUP("Total number of nodes m: %td, n: %td, grid dims L: %g", m, n[0], L[0]);
			// logSETUP("Grid spacing dx: %g, dp: %g", dx, dp);
			// logSETUP("pmin:%g kin_scale: %g", pmin, kin_scale);
			// logSETUP("inv_m:%g inv_nn: %g", inv_m, inv_nn);
		}

		CartesianGrid_(Section& settings)
		{
			double def_dx;
			ind def_n;
			inipp::get_value(settings, "dx", def_dx);
			inipp::get_value(settings, "n", def_n);
			for(DIMS i= 0; i < dims; i++)
			{
				(CoordinateSystem<dims>::dx)[i]= def_dx;
				(CoordinateSystem<dims>::n)[i] = def_n;
			}
			init(n_seq<dims>);
		}

		CartesianGrid_(CoordinateSystem<dims> cs)
			: CoordinateSystem<dims>(cs)
		{
			init(n_seq<dims>);
		}
		CartesianGrid_(double dx, ind n)
		{
			for(DIMS i= 0; i < dims; i++)
			{
				(CoordinateSystem<dims>::dx)[i]= dx;
				(CoordinateSystem<dims>::n)[i] = n;
			}
			init(n_seq<dims>);
			// absorber = AT(this);
		}

		inline cxd scalarProductAll(cxd* from, cxd* to)
		{
			cxd over= 0.;
			for(ind i= 0; i < m; i++) over+= std::conj(to[i]) * from[i];
			return over;
		}

		template<REP R, uind dir> constexpr double pos(ind index)
		{
			if constexpr(R == REP::X) return xmin[dir] + index * dx[dir];
			else
				return dp[dir] *
							 (index < n2[dir] ? index : index - n[dir]);	 // after FFTW 0 freq is at 0 index
		}
	};

	/// @brief Single Cartesian grid divided with MPIDivision across multiple threads
	/// @tparam MPIDivision (default: MPI::Slices)
	/// @tparam dims numbers of dimensions
	template<DIMS dims, class MPIDivision= MPI::Slices>
	using CartesianGrid= CartesianGrid_<dims, MPI::SingleRegion<dims, MPIDivision>>;

	template<DIMS dims, class MPIDivision= MPI::Slices>
	using MultiCartesianGrid= CartesianGrid_<dims, MPI::MultiRegionsReduced<dims, MPIDivision>>;
}		// namespace QSF