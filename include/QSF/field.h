#include "field_p.h"
//probably not needed
struct NoField { static double maxPulseDuration() { return 0.0; } };

template <AXIS Ax, typename ... Pulses>
struct Field : _Operator
{
	static constexpr REP rep = REP::BOTH;
	static constexpr auto axis = Ax;
	using type = Field;
	inline static Field instance;// = nullptr;
	int count = sizeof...(Pulses);
	double delays[sizeof...(Pulses)]{};
	double endtimes[sizeof...(Pulses)]{};
	PulsePrototype* fields[sizeof...(Pulses)]{};

	static constexpr std::string_view name = join_v<Pulses::name...>;

	inline double operator()(double time)
	{
		size_t i = 0;
		double val = 0;
		([&] {
			if (time <= endtimes[i] && time >= delays[i])
				val += static_cast<Pulses*>(fields[i])->operator()(time - delays[i]);
			i++;
		 }(), ...);
		// val = 0.0; i = -1;
		// ((val += ((time <= endtimes[++i] && time >= delays[i]) ? static_cast<Pulses*>(fields[i])->operator()(time - delays[i]) : 0.0)), ...);
		return val;
	}
	double maxPulseDuration()
	{
		size_t i = 0;
		auto max = count > 0 ? endtimes[0] : 0.0;
		for (i = 0; i < count; i++)
			max < endtimes[i] ? max = endtimes[i] : 0.0;
		return max;
	}
	constexpr auto getField(size_t index)
	{
		return static_cast<PulsePrototype*>(fields[index]);
	}

	constexpr Field() : _Operator() {}//= default;
	Field(Pulses&& ... p) :_Operator(),
		delays{ (p.delay_in_cycles * 2.0 * pi / p.omega)... },
		endtimes{ (p.delay_in_cycles * 2.0 * pi / p.omega + p.pulse_time)... },
		fields{ (new Pulses(p))... }
	{
		// logInfo("%p Asssigned to instance", this);
		instance = *this;
	}
	// ~Field() { logInfo("Field Destructor called"); }
	// ~Field() { i = 0; ((delete static_cast<Pulses*>(fields[i++])), ...); }
};

