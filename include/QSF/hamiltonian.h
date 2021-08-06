
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
			if constexpr (R == REP::P) return ((kin_scale[dirs] * coords * coords) + ...);
			else return 0;
		}

		template <REP R, uind ... dirs, typename ... Cooords>
		double potential(Cooords ... coords)
		{
			if constexpr (R == REP::X)
			{
				if constexpr (mcomm.many)
				{
					switch (mcomm.freeCoord)
					{
					case AXIS::NO:
						return _potential.template operator() < AXIS::NO > (coords...);
					case AXIS::X:
						return _potential.template operator() < AXIS::X > (coords...); break;
					case AXIS::Y:
						return _potential.template operator() < AXIS::Y > (coords...); break;
					case AXIS::Z:
						return _potential.template operator() < AXIS::Z > (coords...); break;
					case AXIS::XY:
						return _potential.template operator() < AXIS::XY > (coords...); break;
					case AXIS::YZ:
						return _potential.template operator() < AXIS::YZ > (coords...); break;
					case AXIS::XZ:
						return _potential.template operator() < AXIS::XZ > (coords...); break;
					case AXIS::XYZ:
					default:
						return 0.0; break;
					}
				}
				else return _potential.template operator() < AXIS::NO > (coords...);
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
			else return cos(-val) + I * sin(-val);
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