
namespace QSF
{
/* ::::::::::::::::::::::::::::::: Gauge types :::::::::::::::::::::::::::::::::
Feel free to add more gauges recognized by the physics community. 			  */
	struct VelocityGauge { //p.A
		static constexpr std::string_view CouplingName = "VG";
	};
	struct LengthGauge { // r.F
		static constexpr std::string_view CouplingName = "LG";
	};
	struct CoulombGauge { // Full coupling: p^2/2 - p*A + A^2/2
		static constexpr std::string_view CouplingName = "CG";
	};
	struct MixedGauge { //p.A+r.F
		static constexpr std::string_view CouplingName = "MG";
	};

	/* ::::::::::::::::::::::::::::::: Approx model ::::::::::::::::::::::::::::::::
	Beside gauge type one can use dipole approximation or reduced dimensionality
	models such as models, eg. EckhardtSachaModel which assumes one electron per
	axis in 2D/3D and DipoleApprox(imation). */
	struct NoApprox {
		static constexpr std::string_view CouplingName = "";
	};
	struct DipoleApprox {
		static constexpr std::string_view CouplingName = "DA";
	};
	struct EckhardSachaModel : DipoleApprox {
		static constexpr std::string_view CouplingName = "ES";
	};

	/* :::::::::::::::::::::::::::::::: Couplings ::::::::::::::::::::::::::::::::::
	Important: Gauge types are generic types which do not include any coordinate x/p
	multipliers, as gauges act diffently for different Hamiltonians (if exist).
	Hence the multipliers are applied conditionally at the Hamiltonian level.
	Coupling type does two things:
	1) gathers fields acting on different axesand converts them to potentials
	as needed */
	template <class ApproxModel, class ... Fields>
	struct CouplingBase;

	template <class ... Fields>
	struct CouplingBase<DipoleApprox, Fields...> //: COMPUTATION<double, false, Fields...>
	{
		static_assert(is_unique<Fields...>, "All Fields need to act on different axes!");
		static constexpr uind size = sizeof...(Fields);
		std::tuple<Fields...> fields;
		double last_values[Max(sizeof...(Fields), 1)]{ 0 };
		double last_values_backup[Max(sizeof...(Fields), 1)]{ 0 };

		CouplingBase(Section& settings) : fields(Fields{ settings }...) {}
		CouplingBase(Fields&&...fields) :fields(fields...) {}
		double maxPulseDuration() { return Max(0.0, std::get<Fields>(fields).maxPulseDuration()...); }

		inline void precalc(double time)
		{
			(std::get < Fields>(fields)(time), ...);
		}
		// void reset() { last_values{ 0 }; last_values_backup{ 0 }; }

		// Returns all fields acting on axis ax
		double operator[](AXIS ax)
		{
			// logInfo("operator[]=%g", ((bool(Fields::axis & ax) ? std::get<Fields>(fields).lastVal : 0) + ...));
			return ((bool(Fields::axis & ax) ? std::get<Fields>(fields).lastVal : 0) + ...);
		}

		template <class F> double getValue()
		{
			if constexpr ((std::is_same_v<F, Fields> || ...))
			{
				// logInfo("getValue=%g", std::get<F>(fields).lastVal);
				return std::get<F>(fields).lastVal;
			}
			else return -1; //TODO: log meaningful error 
		}
	};

	template <class ... Fields>
	struct CouplingBase<NoApprox, Fields...>
	{
		std::tuple<Fields...> fields;
		double middle_values[sizeof...(Fields)]{ 0 };
		double middle_values_backup[sizeof...(Fields)]{ 0 };

		static double maxPulseDuration() { return Max(Fields::maxPulseDuration()...); }
		// void reset() { middle_values{ 0 }; middle_values_backup{ 0 }; }

		double operator[](AXIS ax)
		{
			//TODO: NOT IMPLEMENTED!
			return middle_values[ax];
		}
	};

	template <class GaugeType, class ... Fields>
	struct DipoleCoupling;

	//Allows for field auto-promotion to vector potential
	template <class ... Fields>
	struct DipoleCoupling<VelocityGauge, Fields...> :
		CouplingBase<DipoleApprox, std::conditional_t<std::is_base_of_v<_VectorPotential, Fields>, Fields, PromoteFieldToA<Fields>>...>
	{
		using base = CouplingBase<DipoleApprox,
			std::conditional_t<std::is_base_of_v<_VectorPotential, Fields>, Fields, PromoteFieldToA<Fields>>...>;

		static constexpr REP couplesInRep = REP::P;

		DipoleCoupling(Section& settings) : base(settings) {}
		DipoleCoupling(Fields&&...fields) :base(std::forward<Fields>(fields)...) {}

		template <class F> double getValue()
		{
			using target = std::conditional_t<std::is_base_of_v<_VectorPotential, F>, F, PromoteFieldToA<F>>;
			return base::template getValue<target>();
		}
	};

};