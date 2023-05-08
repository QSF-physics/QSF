
#include "field_p.h"

namespace QSF
{
	// probably not needed
	struct NoField
	{
		static double maxPulseDuration() { return 0.0; }
	};

	template<AXIS Ax, typename... Pulses> struct Field: _Operator
	{
		using type								= Field;
		static constexpr REP rep	= REP::NONE;
		static constexpr bool late= false;

		static constexpr auto axis						= Ax;
		static constexpr std::string_view name= join_v<Pulses::name...>;
		static constexpr uind count						= sizeof...(Pulses);
		std::tuple<Pulses...> pulses;
		double delays[sizeof...(Pulses)]{};
		double endtimes[sizeof...(Pulses)]{};
		double lastVal= 0;
		inline double operator()(double time)
		{
			size_t i= 0;
			lastVal = 0;
			std::apply(
					[&](auto const&... pulse)
					{
						(((lastVal+=
							 (!(time > endtimes[i] || time < delays[i]) ? pulse(time - delays[i]) : 0.0)),
							i++),
						 ...);
					},
					pulses);
			return lastVal;
		}
		double maxPulseDuration()
		{
			size_t i= 0;
			auto max= count > 0 ? endtimes[0] : 0.0;
			for(i= 0; i < count; i++) max < endtimes[i] ? max= endtimes[i] : 0.0;
			return max;
		}

		Field(Pulses&&... p)
				: pulses{std::forward<Pulses>(p)...},
					delays{(p.delay_in_cycles * 2.0 * pi / p.omega)...},
					endtimes{(p.delay_in_cycles * 2.0 * pi / p.omega + p.pulse_time)...}
		{}
		Field(Section& settings)
				: pulses{Pulses{PulsePrototype{settings, axisName(Ax)}}...},
					delays{(std::get<Pulses>(pulses).delay_in_cycles * 2.0 * pi /
									std::get<Pulses>(pulses).omega)...},
					endtimes{(std::get<Pulses>(pulses).delay_in_cycles * 2.0 * pi /
												std::get<Pulses>(pulses).omega +
										std::get<Pulses>(pulses).pulse_time)...}
		{}
	};

	struct _VectorPotential
	{};

	template<AXIS Ax, typename... Pulses> struct VectorPotential: _Operator, _VectorPotential
	{
		using type								= VectorPotential;
		static constexpr REP rep	= REP::NONE;
		static constexpr bool late= false;

		static constexpr auto axis						= Ax;
		static constexpr std::string_view name= join_v<Pulses::name...>;
		static constexpr uind count						= sizeof...(Pulses);
		std::tuple<Pulses...> pulses;
		double delays[sizeof...(Pulses)]{};
		double endtimes[sizeof...(Pulses)]{};
		double lastVal= 0;
		inline double operator()(double time)
		{
			size_t i= 0;
			lastVal = 0;
			std::apply(
					[&](auto const&... pulse)
					{
						(((lastVal+=
							 (!(time > endtimes[i] || time < delays[i]) ? pulse(time - delays[i]) : 0.0)),
							i++),
						 ...);
					},
					pulses);
			return lastVal;
		}
		double maxPulseDuration()
		{
			size_t i= 0;
			auto max= count > 0 ? endtimes[0] : 0.0;
			for(i= 0; i < count; i++) max < endtimes[i] ? max= endtimes[i] : 0.0;
			return max;
		}

		VectorPotential(Pulses&&... p)
				: pulses{std::forward<Pulses>(p)...},
					delays{(p.delay_in_cycles * 2.0 * pi / p.omega)...},
					endtimes{(p.delay_in_cycles * 2.0 * pi / p.omega + p.pulse_time)...}
		{}
		VectorPotential(Section& settings)
				: pulses{Pulses{PulsePrototype{settings, axisName(Ax)}}...},
					delays{(std::get<Pulses>(pulses).delay_in_cycles * 2.0 * pi /
									std::get<Pulses>(pulses).omega)...},
					endtimes{(std::get<Pulses>(pulses).delay_in_cycles * 2.0 * pi /
												std::get<Pulses>(pulses).omega +
										std::get<Pulses>(pulses).pulse_time)...}
		{}
	};

	template<class F> struct PromoteFieldToA: F, _VectorPotential
	{
		using type= PromoteFieldToA;
		PromoteFieldToA(F&& field): F{std::forward<F>(field)} {}
		PromoteFieldToA(Section& settings): F(settings) {}
		using F::lastVal;

		// double Fval = 0.0;
		double prevVal = 0.0;
		double prevTime= 0.0;
		inline double operator()(double time)
		{
			prevVal= lastVal;
			F::operator()(time);

			// logInfo("precalc: F(t)=%g A(t)=%g ", lastVal, prevVal + (prevTime - time) * lastVal);
			lastVal=
					prevVal + (prevTime - time) * lastVal;	 // includes the minus in A(t)=int_0^t -F(t) dt
			prevTime= time;
			return lastVal;
		}
		// static constexpr REP rep = F::rep;
		// static constexpr bool late = false;
		// static constexpr auto axis = F::Ax;
		// static constexpr std::string_view name = F::name;
		// static constexpr uind count = F::count;

		// double delays[sizeof...(Pulses)]{};
		// double endtimes[sizeof...(Pulses)]{};
		// double lastVal = 0;
	};
}		// namespace QSF