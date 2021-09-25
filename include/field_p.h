namespace QSF
{


// DATATYPE FOR DESIGNATED INITIALIZATION
	struct _FieldConfig
	{
		/* If you think that a generic field parameter is lacking add it here */
		double field, omega, ncycles, FWHM_percent, phase_in_pi_units, delay_in_cycles;
		double Up() const { return 0.25 * field / omega / omega; }
	};

	//Sets a few derived quantities for convienience 
	struct PulsePrototype : _FieldConfig
	{
		double phi, pulse_time, onecycle_time;
		// PulsePrototype() = default;
		PulsePrototype(Section& settings, std::string dir = "") : _FieldConfig()
		{
			dir = "::" + dir;
			inipp::get_value(settings, "F" + dir, field);
			inipp::get_value(settings, "omega" + dir, omega);
			inipp::get_value(settings, "ncycles" + dir, ncycles);
			inipp::get_value(settings, "FWHM" + dir, FWHM_percent);
			inipp::get_value(settings, "phase" + dir, phase_in_pi_units);
			inipp::get_value(settings, "delay" + dir, delay_in_cycles);
		}

		PulsePrototype(_FieldConfig p) : _FieldConfig(p),
			phi(p.phase_in_pi_units* pi),
			pulse_time(2.0 * pi * p.ncycles / p.omega),
			onecycle_time(2.0 * pi / p.omega) {}
	};

	// PREDEFINED PULSES
	//TODO: Implement extra step after ntsteps fixing vector potential to exact 0

	////////////////////////////////////////
	// PREDEFINED FIELDS ///////////////////
	struct EmptyPulse :PulsePrototype
	{
		static constexpr std::string_view name = "EMPTY";
		// EmptyPulse() = default;
		EmptyPulse(_FieldConfig p) : PulsePrototype(p) {}
	};
	struct SinPulse :PulsePrototype
	{
		static constexpr std::string_view name = "SIN";
		// SinPulse() = default;
		SinPulse(_FieldConfig p) :PulsePrototype(p) {}
		inline double operator()(double time) const
		{
			return field * sin(omega * time + phi);
		}
	};
	struct ConstantPulse :PulsePrototype
	{
		static constexpr std::string_view name = "--";
		// ConstantPulse() = default;
		ConstantPulse(_FieldConfig p) :PulsePrototype(p) {}
		inline double operator()(double time) const
		{
			return field;
		}
	};
	struct ChemPhysPulse :PulsePrototype
	{
		static constexpr std::string_view name = "CHEM";
		double env_mult;
		// ChemPhysPulse() = default;
		ChemPhysPulse(_FieldConfig p) :PulsePrototype(p),
			env_mult(p.omega / (2.0 * p.ncycles)) {}
		inline double operator()(double time) const
		{
			double tmp = omega * time + phi - pi * ncycles;
			return field * (sin(env_mult * time) * cos(tmp) + cos(env_mult * time) * sin(tmp) / ncycles);
		}
	};

	//PREDEFINED ENVELOPES

	//should be used together with ChemPhysPulse
	template <typename Pulse>
	struct ChemPhysEnvelope : Pulse
	{
		static constexpr std::string_view _name = "ENV_";
		static constexpr std::string_view name = join_v<_name, Pulse::name>;
		double env_mult;
		// ChemPhysEnvelope() = default;
		ChemPhysEnvelope(_FieldConfig p) : Pulse{ p }, env_mult(p.omega / (2.0 * p.ncycles)){}
		inline double operator()(double time) const
		{
			if constexpr (std::is_same_v<Pulse, EmptyPulse>) return sin(env_mult * time);
			else return  Pulse::operator()(time) * sin(env_mult * time);
		}
	};

	//Ignores FWHM
	template <typename Pulse>
	struct Sin2Envelope : Pulse
	{

		static constexpr std::string_view _name = "SIN2_";
		static constexpr std::string_view name = join_v<_name, Pulse::name>;
		double env_mult;
		// Sin2Envelope() = default;
		Sin2Envelope(_FieldConfig p) : Pulse{ p }, env_mult(p.omega / (2.0 * p.ncycles)) {
			//TODO: Warn if FWHM >0
		}
		inline double operator()(double time) const
		{
			if constexpr (std::is_same_v<Pulse, EmptyPulse>) return Power(sin(env_mult * time), 2);
			else return Pulse::operator()(time) * Power(sin(env_mult * time), 2);
		}
	};

	template <typename Pulse>
	struct RampEnvelope : Pulse
	{
		static constexpr std::string_view _name = "RAMP_";
		static constexpr std::string_view name = join_v<_name, Pulse::name>;
		double leftEdge;
		double rightEdge;
		double ramp;
		// RampEnvelope() = default;
		RampEnvelope(_FieldConfig p) : Pulse{ p },
			leftEdge((-2.0 * (-1.0 + p.FWHM_percent))* Pulse::pulse_time),
			rightEdge((-1.0 + 2.0 * p.FWHM_percent)* Pulse::pulse_time),
			ramp(0.5 * (1.0 / (Pulse::pulse_time * (1.0 - p.FWHM_percent))))
		{
			//TODO: Throw. supported FWHM range: [0.75-1]
		};
		inline double operator()(double time) const
		{
			if constexpr (std::is_same_v<Pulse, EmptyPulse>)
			{
				if (time < leftEdge) return ramp * time;
				else if (time > rightEdge) return ramp * (Pulse::pulse_time - time);
				else return 1.0;
			}
			else
			{
				if (time < leftEdge) return Pulse::operator()(time) * ramp * time;
				else if (time > rightEdge) return Pulse::operator()(time) * ramp * (Pulse::pulse_time - time);
				else return Pulse::operator()(time) * 1.0;
			}
		}
	};

	/* This envelope sets the power of sine
	such that FWHM_percent is satisfied */
	template <typename Pulse>
	struct SinEnvelope : Pulse
	{
		static constexpr std::string_view _name = "ENV_";
		static constexpr std::string_view name = join_v<_name, Pulse::name>;
		double sin_pow;
		double env_mult;
		// SinEnvelope() = default;
		SinEnvelope(_FieldConfig p) : Pulse{ p },
			env_mult(p.omega / (2.0 * p.ncycles)),
			sin_pow(log(2) / log(1 / sin(0.5 * pi * (1 - p.FWHM_percent)))){
				// Supported FWHM range: (0,1)
		}
		inline double operator()(double time) const
		{
			if constexpr (std::is_same_v<Pulse, EmptyPulse>)
				return pow(sin(env_mult * time), sin_pow);
			else
				return Pulse::operator()(time) * pow(sin(env_mult * time), sin_pow);
		}
	};
	template <typename Pulse>
	struct GaussianEnvelope : Pulse
	{
		static constexpr std::string_view _name = "GAUSS_";
		static constexpr std::string_view name = join_v<_name, Pulse::name>;
		double gauss_pow;
		// GaussianEnvelope() = default;
		GaussianEnvelope(_FieldConfig p) :
			Pulse{ p },
			gauss_pow(-4.0 * log(2) / Power(p.FWHM_percent * Pulse::pulse_time, 2))
		{
			//TODO: Warn. Supported FWHM range: (0.0-0.35], otherwise might not start from zero
		}
		inline double operator()(double time) const
		{
			if constexpr (std::is_same_v<Pulse, EmptyPulse>)
				return exp(gauss_pow * Power(time - 0.5 * Pulse::pulse_time, 2));
			else
				return Pulse::operator()(time) * exp(gauss_pow * Power(time - 0.5 * Pulse::pulse_time, 2));
		}
	};
};