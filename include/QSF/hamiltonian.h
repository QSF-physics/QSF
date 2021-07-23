
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
		using Base::n0_l;
		using Base::n0_o;
		using Base::DIM;
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
		template <REP R, typename ... Args>
		double couple(Args ... args)
		{
			return 0.0;
			// if constexpr (C_Op::couplesInRep == R && C_Op::size)
			// 	return ((args * coupling[AXIS::NO]) + ...);
			// else return 0.0;
		}
		template < REP R, typename ... Cooords>
		inline double operator()(Cooords ... coords)
		{
			if constexpr (R == REP::P)
				return (InducedGrid::kin_scale * ((coords * coords) + ...)) - couple<R>(coords...);
			else return potential(coords...) + couple<R>(coords...);
		}

		template <class Op, class...Coords>
		inline double call(Coords...coords)
		{
			// logInfo("call");
			if constexpr (std::is_same_v<Op, KineticEnergy>)
				return operator() < REP::P > (coords...);
			else if constexpr (std::is_same_v<Op, PotentialEnergy>)
				return operator() < REP::X > (coords...);
			else return 3.0;
		}


		template <MODE M, REP R>
		auto expOp(double val)
		{
			if constexpr (M == MODE::IM)
				return exp(-val);
			else
				return cos(-val) + I * sin(-val);
		}

	};
}




// template < REP R>
// inline double operator()(double x, double y)
// {
// 	if constexpr (R == REP::P)
// 		return (InducedGrid::kin_scale * (x * x * x * x + y * y));
// 	else return potential(x * x * x * x + y * y);
// }