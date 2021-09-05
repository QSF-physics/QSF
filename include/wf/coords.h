#pragma once
namespace QSF
{
	template <DIMS size>
	struct CoordinateSystem {
		double dx[size];
		ind n[size];
		static constexpr DIMS DIM = size;
	};


	template <DIMS size, class MPIRegions>
	struct CartesianGrid_ :CoordinateSystem<size>
	{
		// static_assert(std::is_base_of_v< AbsorberType, AT>,
		// "Second CartesianGrid template argument must be type derived from AbsorberType");

		using MPIGridComm = MPIRegions;
		using MPIDivision = typename MPIRegions::MPIDivision;
		/* Init of the above part is postponed to the wf */

		using CoordinateSystem<size>::n;
		using CoordinateSystem<size>::dx;
		using CoordinateSystem<size>::DIM;

		ind n2[size];
		ind nn[size];
		ind m;
		// double inv_nn;

		double inv_m;
		double dV;
		double dVm;
		double L[size];
		double xmin[size];
		double inv_dx[size];
		double inv_2dx[size];
		double inv_n[size];

		double dp[size];
		double pmin[size];
		double pmax[size];
		double kin_scale[size];
		double dVP;

		template <REP R>
		double vol()
		{
			if constexpr (R == REP::X) return dV;
			else return dVm;
		}

		template <uind ...Is>
		void init(seq<Is...>)
		{
			m = 1;
			([&] {
				n2[Is] = n[Is] / 2;
				nn[Is] = n[Is] * n[Is];
				inv_n[Is] = 1.0 / n[Is];
				m *= n[Is];
			 }(), ...);
			inv_m = 1.0 / m;

			dV = 1.0;
			dVP = 1.0;
			([&] {
				dV *= dx[Is];
				L[Is] = (dx[Is] * (n[Is] - 1));
				xmin[Is] = -0.5 * L[Is];
				inv_dx[Is] = 1.0 / dx[Is];
				inv_2dx[Is] = inv_dx[Is] / 2.0;
				dp[Is] = 2.0 * pi / double(n[Is]) / dx[Is];
				pmin[Is] = -pi / dx[Is];
				pmax[Is] = -pmin[Is] - dp[Is];
				kin_scale[Is] = dp[Is] * dp[Is] * 0.5;
				dVP *= dp[Is];
			 }(), ...);
			dVm = dV * inv_m;

			// logSETUP("CartesianGrid init");
			// logSETUP("Total number of nodes m: %td, n: %td, grid size L: %g", m, m, L);
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
			for (DIMS i = 0; i < size; i++)
			{
				(CoordinateSystem<size>::dx)[i] = def_dx;
				(CoordinateSystem<size>::n)[i] = def_n;
			}
			init(n_seq<size>);

		}

		CartesianGrid_(CoordinateSystem<size> cs) : CoordinateSystem<size>(cs)
		{
			init(n_seq<size>);
		}
		CartesianGrid_(double dx, ind n)
		{
			for (DIMS i = 0; i < size; i++)
			{
				(CoordinateSystem<size>::dx)[i] = dx;
				(CoordinateSystem<size>::n)[i] = n;
			}
			init(n_seq<size>);
			// absorber = AT(this);
		}

		inline cxd scalarProductAll(cxd* from, cxd* to)
		{
			cxd over = 0.;
			for (ind i = 0; i < m; i++)
				over += std::conj(to[i]) * from[i];
			return over;
		}

		template <REP R, uind dir>
		constexpr double pos(ind index)
		{
			if constexpr (R == REP::X) return xmin[dir] + index * dx[dir];
			else return dp[dir] * (index >= n2[dir] ? index - n[dir] : index); //after FFTW 0 freq is at 0 index
		}
	};


	template <DIMS size, class MPIDivision = MPI::Slices>
	using CartesianGrid = CartesianGrid_<size, MPI::SingleRegion<size, MPIDivision>>;

	template <DIMS size, class MPIDivision = MPI::Slices>
	using MultiCartesianGrid = CartesianGrid_<size, MPI::MultiRegionsReduced<size, MPIDivision>>;
}