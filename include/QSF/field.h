#include "field_p.h"
//probably not needed
struct NoField { static double maxPulseDuration() { return 0.0; } };

template <AXIS Ax, typename ... Pulses>
struct Field : _Operator
{
	using type = Field;
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = false;

	static constexpr auto axis = Ax;
	static constexpr std::string_view name = join_v<Pulses::name...>;
	static constexpr uind count = sizeof...(Pulses);
	std::tuple<Pulses...> pulses;
	double delays[sizeof...(Pulses)]{};
	double endtimes[sizeof...(Pulses)]{};
	double lastVal = 0;
	inline double operator()(double time)
	{
		size_t i = 0;
		lastVal = 0;
		std::apply
		(
			[&](auto const&... pulse)
			{
				((lastVal += (!(time > endtimes[i] || time < delays[i++]) ?
							  pulse(time - delays[i]) : 0.0)), ...);
							  // if (time <= endtimes[i] && time >= delays[i])
							  // 	val += field.operator()(time - delays[i]);
							  // i++;
			}, pulses
		);

		// ([&] {
		// 	if (time <= endtimes[i] && time >= delays[i])
		// 		val += static_cast<Pulses*>(pulses[i])->operator()(time - delays[i]);
		// 	i++;
		//  }(), ...);
		// val = 0.0; i = -1;
		// ((val += ((time <= endtimes[++i] && time >= delays[i]) ? static_cast<Pulses*>(pulses[i])->operator()(time - delays[i]) : 0.0)), ...);
		return lastVal;
	}
	double maxPulseDuration()
	{
		size_t i = 0;
		auto max = count > 0 ? endtimes[0] : 0.0;
		for (i = 0; i < count; i++)
			max < endtimes[i] ? max = endtimes[i] : 0.0;
		return max;
	}
	// constexpr auto getField(size_t index)
	// {
	// 	return static_cast<PulsePrototype*>(pulses[index]);
	// }

	Field(Pulses&& ... p) :
		pulses{ std::move(p)... },
		delays{ (p.delay_in_cycles * 2.0 * pi / p.omega)... },
		endtimes{ (p.delay_in_cycles * 2.0 * pi / p.omega + p.pulse_time)... }{}
	// ~Field() { logInfo("Field Destructor called"); }
	// ~Field() { i = 0; ((delete static_cast<Pulses*>(fields[i++])), ...); }
};

