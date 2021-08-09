
namespace Schrodinger
{
	//CRTP 
	template <class GType, class V_Op, class C_Op = DipoleCoupling<VelocityGauge> >
	struct Spin0 : LocalGrid < Spin0<GType, V_Op, C_Op>, GType, 1>
	{
		using Base = LocalGrid< Spin0<GType, V_Op, C_Op>, GType, 1>;
		using Base::psi;
		using Base::postCompute;
		using Base::DIM;
		using Base::kin_scale;
		using Base::mcomm;
		using Base::MPIGridComm;
		static constexpr REP couplesInRep = C_Op::couplesInRep;

		V_Op _potential;
		C_Op _coupling;

		Spin0(Section& settings) : Base(settings), _potential(settings), _coupling(settings)
		{
			logInfo("Spin0 init");
		}

		Spin0(GType gtype, V_Op potential, C_Op coupling) :Base(gtype), _potential(potential), _coupling(coupling) {}

		// static constexpr REP rep = REP::BOTH;
		// static constexpr string_view name = "Schrodinger";

		template <REP R, uind ... dirs, typename ... Cooords>
		double coupling(Cooords ... coords)
		{
			if constexpr (C_Op::couplesInRep == R && C_Op::size)
				return ((coords * _coupling[Axis<dirs>]) + ...);
			else return 0.0;
		}

		template <REP R, uind ... dirs, typename ... Cooords>
		double kinetic(Cooords ... coords)
		{
			if constexpr (R == REP::P)
				return ((0.5 * coords * coords) + ...);
				// return ((kin_scale[dirs] * coords * coords) + ...);
			else return 0;
		}

		template <REP R, uind ... dirs, typename ... Cooords>
		double potential(Cooords ... coords)
		{
			if constexpr (R == REP::X)
			{
				if constexpr (Base::MPIGridComm::many)
					return _potential.partial(mcomm.freeCoord, coords...);
				else return _potential.operator()(coords...);
			}
			else return 0;
		}

		template < REP R, uind ... dirs, typename ... Cooords>
		inline double operator()(Cooords ... coords)
		{
			if constexpr (R == REP::P)
				return kinetic<R, dirs...>(coords...) - coupling<R, dirs...>(coords...);
			else return potential<R, dirs...>(coords...) + coupling<R, dirs...>(coords...);
		}

		template <class Op, uind ... dirs, class...Coords>
		inline double call(Coords...coords)
		{
			if constexpr (std::is_same_v<Op, Identity>)
				return 1.0;
			else if constexpr (std::is_same_v<Op, KineticEnergy>)
				return kinetic<REP::P, dirs...>(coords...);
			else if constexpr (std::is_same_v<Op, PotentialEnergy>)
				return potential<REP::X, dirs...>(coords...);
			// else 
			return 3.0;
		}


		template <MODE M, REP R>
		auto expOp(double val)
		{
			if constexpr (M == MODE::IM) return exp(-val);
			else return cxd{ cos(-val) , sin(-val) };
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