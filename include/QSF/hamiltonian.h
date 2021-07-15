struct KineticEnergy {
	static constexpr REP rep = REP::P;
};
struct PotentialEnergy {
	static constexpr REP rep = REP::X;
};
struct Identity {
	static constexpr REP rep = REP::NONE;
};
struct TotalEnergy
{
	static constexpr REP rep = REP::NONE;
};
struct EnergyDifference
{
	static constexpr REP rep = REP::NONE;
};
struct Symmetrize
{
	static constexpr REP rep = REP::NONE;
};
struct AntiSymmetrize
{
	static constexpr REP rep = REP::NONE;
};
struct Orthogonalize
{
	static constexpr REP rep = REP::NONE;
};
struct Normalize
{
	static constexpr REP rep = REP::NONE;
};

namespace Schrodinger
{
	/* Define named operators that can be used in computations */

	struct Couplings { };

	template <class V_Op, class C_Op, class GType>
	struct Spin0 : V_Op, C_Op, WF < Spin0<V_Op, C_Op, GType>, GType, 1>
	{
		using wf = WF < Spin0<V_Op, C_Op, GType>, GType, 1>;
		using wf::psi;
		using wf::post_evolve;
		using grid = typename wf::grid;
		// using wf::transfer;
		using V = V_Op;
		using C = C_Op;

		Spin0(Section& settings) : wf(settings)
		{
			logInfo("Spin0 init");
		}
		static constexpr REP couplesInRep = C_Op::couplesInRep;
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


		template <AXIS AX, class Op>
		inline double avg()
		{
			return 1.0;
		}
		// template <class T>
		// auto compute();

		// template <class Op>
		// auto compute(CO_AVG<Op>)
		// {
		// 	static_assert(1, "operation not implemented");
		// }

		// template <AXIS AX>
		// auto compute(AVG<AX, KinEnergy>)
		// {
		// 	return operator() < REP::P, OPTIMS::NONE > (0, 0, 0);
		// }
		// template <AXIS AX>
		// auto compute(AVG<AX, PotEnergy>)
		// {
		// 	return operator() < REP::P, OPTIMS::NONE > (0, 0, 0);
		// }

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