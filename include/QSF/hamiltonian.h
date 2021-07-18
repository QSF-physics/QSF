
namespace Schrodinger
{
	/* Define named operators that can be used in computations */

	struct Couplings { };

	template <class GType, class V_Op, class C_Op = DipoleCoupling<VelocityGauge> >
	struct Spin0 : WF < Spin0<GType, V_Op, C_Op>, GType, 1>
	{
		using Base = WF < Spin0<GType, V_Op, C_Op>, GType, 1>;
		using Base::psi;
		using Base::post_step;
		using InducedGrid = typename Base::InducedGrid;
		static constexpr REP couplesInRep = C_Op::couplesInRep;

		V_Op potential;
		C_Op coupling;

		Spin0(Section& settings) : Base(settings), potential(settings), coupling(settings)
		{
			logInfo("Spin0 init");
		}

		Spin0(GType gtype, V_Op potential, C_Op coupling) :Base(gtype), potential(potential), coupling(coupling) {}

		// V pot;
		// static constexpr REP rep = REP::BOTH;
		// static constexpr string_view name = "Schrodinger";

		template < REP R, typename ... Args>
		inline double operator()(Args ... args)
		{
			if constexpr (R == REP::P)
			{
				if constexpr (C_Op::couplesInRep == R && C_Op::size)
					return (InducedGrid::kin_scale * (Power(args, 2) + ...) +
							C_Op::template operator() < R, OPTIMS::NONE > (args...));
				else return InducedGrid::kin_scale * (Power(args, 2) + ...);
			}
			else if (R == REP::X)
			{
				if constexpr (C_Op::couplesInRep == R && C_Op::size)
					return (potential(args...) +
							C_Op::template operator() < R, OPTIMS::NONE > (args...));
				else return potential(args...);
			}
		}


		template <class Op, class...Coords>
		inline double call(Coords...coords)
		{
			// logInfo("call");
			if constexpr (std::is_same_v<Op, KineticEnergy>)
				return operator() < REP::P > (coords...);
			else if constexpr (std::is_same_v<Op, PotentialEnergy>)
			{
				return 1.0;//operator() < REP::X > (coords...);
			}
			else //if constexpr (std::is_same_v<Op, KineticEnergy>)
				return 3.0;
		}


		template <MODE M>
		auto expOp(double val)
		{
			if constexpr (M == MODE::IM)
				return exp(-val);
			else
				return cos(-val) + I * sin(-val);
		}
		using Base::local_n;
		using Base::local_start;
		using Base::DIM;

		template <MODE M, REP R>
		void evolve(double delta)
		{
			// for (ind i = 0; i < local_n; i++)
			// {
			// 	if constexpr (DIM == 1)
			// 		psi[i] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (i + local_start));
			// 	else
			// 	{
			// 		ind readInd1 = i * InducedGrid::n;
			// 		//Due to FFTW flag FFTW_MPI_TRANSPOSED_OUT we need to switch x<->y for DIM>1
			// 		for (ind j = 0; j < InducedGrid::n; j++)
			// 		{
			// 			ind readInd2 = readInd1 + j;
			// 			if constexpr (DIM == 2)
			// 			{
			// 				if (R == REP::X || MPI::region)
			// 					psi[readInd2] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (i + local_start, j));
			// 				else
			// 					psi[readInd2] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (j, i + local_start));
			// 			}
			// 			else
			// 			{
			// 				ind readInd2 = readInd2 * InducedGrid::n;
			// 				for (ind k = 0; k < InducedGrid::n; k++)
			// 				{
			// 					ind readInd3 = readInd2 + k;
			// 					if (R == REP::X || MPI::region)
			// 						psi[readInd3] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (i + local_start, j, k));
			// 					else
			// 						psi[readInd3] *= expOp<M>(delta * operator() < R, OPTIMS::NONE > (j, i + local_start, k));
			// 				}
			// 			}
			// 		}
			// 	}
			// }
		}

	};
}