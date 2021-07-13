
namespace Schrodinger
{
	/* Define named operators that can be used in computations */
	struct KinEnergy { };
	struct PotEnergy { };
	struct Couplings { };

	template <class V_Op, class C_Op, class GType>
	struct Spin0 : V_Op, C_Op, WF < Spin0<V_Op, C_Op, GType>, GType, 1>
	{
		using wf = WF < Spin0<V_Op, C_Op, GType>, GType, 1>;
		using wf::psi;
		using grid = typename wf::grid;
		using V = V_Op;
		using C = C_Op;
		Spin0() {}
		static constexpr REP couplesInRep = C::couplesInRep;
		// V pot;
		// static constexpr REP rep = REP::BOTH;
		// static constexpr string_view name = "Schrodinger";

		template < REP R, OPTIMS opt, typename ... Args>
		double operator()(Args ... args) const
		{
			if constexpr (R == REP::P)
			{
				if constexpr (couplesInRep == R)
					return (grid::kin_scale * (Power(grid::pos<R>(args), 2) + ...) +
							C_Op::template operator() < R, OPTIMS::NONE > (args...));
				else return grid::kin_scale * (Power(grid::pos<R>(args), 2) + ...);
			}
			else if (R == REP::X)
			{
				if constexpr (couplesInRep == R)
					return (V::template operator() < R, opt > (grid::pos<R>(args)...) +
							C_Op::template operator() < R, OPTIMS::NONE > (grid::pos<R>(args)...));
				else return V::template operator() < R, opt > (grid::pos<R>(args)...);
			}
		}
		template <MODE M>
		auto expOp(double val)
		{
			if constexpr (M == MODE::IM)
				return exp(-val);
			else
				return cos(-val) + I * sin(-val);
		}
		using wf::local_n;
		using wf::local_start;
		using wf::DIM;

		template <MODE M, REP R>
		void evolve(double delta)
		{
			// for (ind i = 0; i < local_n; i++)
			// {
			// 	if constexpr (DIM == 1)
			// 		psi[i] *= expOp(delta * operator() < R, OPTIMS::NONE > (i + local_start));
			// 	else
			// 	{
			// 		readInd1 = i * grid::n;
			// 		//Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to switch x<->y for DIM>1
			// 		for (j = 0; j < grid::n; j++)
			// 		{
			// 			readInd2 = readInd1 + j;
			// 			if constexpr (DIM == 2)
			// 			{
			// 				if (R == REP::X || MPI::region)
			// 					psi[readInd2] *= te(delta * operator() < R, OPTIMS::NONE > (i + local_start, j));
			// 				else
			// 					psi[readInd2] *= te(delta * operator() < R, OPTIMS::NONE > (j, i + local_start));
			// 			}
			// 			else
			// 			{
			// 				readInd2 = readInd2 * grid::n;
			// 				for (k = 0; k < grid::n; k++)
			// 				{
			// 					readInd3 = readInd2 + k;
			// 					if (R == REP::X || MPI::region)
			// 						psi[readInd3] *= te(delta * operator() < R, OPTIMS::NONE > (i + local_start, j, k));
			// 					else
			// 						psi[readInd3] *= te(delta * operator() < R, OPTIMS::NONE > (j, i + local_start, k));
			// 				}
			// 			}
			// 		}
			// 	}
			// }
		}
	};
}