template <class PASSES, class PROP, class Worker>
struct Routine;

template <uind ...PASS, class PROP, class Worker>
struct Routine<seq<PASS...>, PROP, Worker>
{
	static constexpr auto mode = PROP::mode;

	PROP propagator;

	Worker worker;
	std::string_view name;


	explicit Routine(PROP propagator) : propagator(propagator)
	{
	};

	Routine(Worker&& worker) :
		worker(std::move(worker))
	{
		name = PROP::name;

	}
	// LambdaClass(const F &lambdaFunc_): lambdaFunc(lambdaFunc_) {}
	// LambdaClass(F &&lambdaFunc_) : lambdaFunc(std::move(lambdaFunc_)) {}
	// LambdaClass(F& _lambdaFunc) : lambdaFunc(_lambdaFunc) {}
	// LambdaClass(F&& _lambdaFunc) : lambdaFunc(std::forward<F>(_lambdaFunc)) {}
	void greet()
	{
		logImportant("EXECUTING ROUTINE [%s] IN MODE [%s] REGIONS: [%d] USING OPERATOR SPLIT: [%s]",
					 name.data(),
					 modeName(mode),
					 MPI::regionCount,
					 PROP::name.data());
	}


	void setupOptimizations()
	{
	// 	constexpr auto rout = get<RI>(ROUTINES);
	// 	if (andQ<RI>(PRECOMP_COORDS)) precomputeCoords();
	// 	initPulseOptimizations<RI>();
	// 	if (andQ<RI>(C_DIPACC_X)) buildAndScatter(makeVStatDer<AXIS::X>, opt_dvstat_dx);
	// 	if (andQ<RI>(C_DIPACC_Y) && DIM > 1) buildAndScatter(makeVStatDer<AXIS::Y>, opt_dvstat_dy);
	// 	if (andQ<RI>(C_DIPACC_Z) && DIM > 2) buildAndScatter(makeVStatDer<AXIS::Z>, opt_dvstat_dz);
	}

	void run()
	{
		greet();
		/* configChecks<RI>();
		config();
		MPI::test();
		propagator.setup();
		wf.setup();
		Eigen::setup(eigenUsedByComp<typename RT::template Output_t<0>>, imaginaryTimeQ);
		constexpr auto seq = (typename RT::Passes_t){}; */
		(run(PASS), ...);
	}

	~Routine()
	{
		/* logInfo("RoutineEnv destructs");
		Eigen::destroy();
		// WAVE::destroySlice();
		if (psi_total != nullptr)
		{
			delete[] psi_total;
			psi_total = nullptr;
		}
		// // closeFile(file_log);
		// Eigen::destroy();

		// fftw_free(WF::psi);
		// fftw_destroy_plan(WF::transf_p2x);
		// fftw_destroy_plan(transf_p2x);
		// http://www.fftw.org/fftw3_doc/MPI-Initialization.html
		fftw_mpi_cleanup(); */
	}

	template <typename WHEN>
	inline void reapData()
	{
		// runDumps<RI, R, WHEN>();
		// Dumps_t::template run<startREP, WHEN>();

		// 	propagator.template fourier<REP::BOTH^ startREP>();
		// 	Dumps_t::template run< REP::BOTH^ startREP, WHEN>();
		// 	outputs.template run < mode, LATE, REP::BOTH^ startREP, optims>(propagator);
		// 	propagator.template fourier<startREP>();
		// 	outputs.template logOrPass<WHEN>(1);
	}

	inline void start()
	{

	}

	void run(uind i)
	{
		propagator.run(std::move(worker), i);
		// // Timings::start(name, RI, PASS);
		// worker(ind(WHEN::AT_START), i, propagator.wf);
		// // if (eigenUsedByCompQ<Output_t>)
		// // 	Eigen::load<RI, PASS, mode>(eigenUsedByComp<Output_t>);
		// propagator.reset();
		// if (PROP::calcsEnabled())
		// {
		// 	// setupInitialState<PASS, REP::X, mode>();
		// 	// outputs.setupComputations();
		// }
		// // if constexpr (startREP == REP::P)
		// // {
		// // 	logInfo("Starts in P");
		// // 	propagator.template fourier<REP::P>();
		// // }
		// propagator.outputs.writeCaptions();

		// while (propagator.stillEvolving());
		// {
		// 	propagator.makeStep((typename PROP::ChainExpander) {});
		// 	// reapData<DURING<>>();
		// 	worker(ind(propagator.step), i, propagator.wf);
		// }
		// worker(ind(WHEN::AT_END), i, propagator.wf);

		// Timings::stop();
	}

};
template <class PASSES, class PROP, class Worker>
Routine<PASSES, PROP, Worker> Repeat(Worker&& w)
{
	return Routine<PASSES, PROP, Worker>(std::move(w));
}


// im load
/*
void loadIMConfig()
{
	using namespace inipp;
	get_value(section, "n", n);
	get_value(section, "L", L);
	get_value(section, "max_imaginary_steps", max_imaginary_steps);
	get_value(section, "state_accuracy", state_accuracy);

	dx = L / double(n - 1);
	dV = pow(dx, DIM);
	sqrt_dV = sqrt(dV);
	inv_dx = 1.0 / dx;
	inv_2dx = inv_dx / 2.0;
	Im::n = n;
	Im::L = L;
	Im::nn = POW2(Im::n);
	Im::m = pow(Im::n, DIM);
}

template <ind SRC>
void loadSourceConfig()
{
	using namespace inipp;
	logImportant("SOURCE ROUTINE INDEX: [%zu]", SRC);
	auto sectionName = "";//modeName(Routine_t<SRC>::mode);
	auto section = ini.sections[sectionName];
	get_value(section, "n", Src::n);
	//Use the overload with name if present
	auto sectionName2 = "";//Routine_t<SRC>::name;
	section = ini.sections[sectionName2];
	get_value(section, "n", Src::n);
}
void loadFromSection(string_view nameR)
{
	if (nameR.length() == 0) return;
	using namespace inipp;
	auto section = ini.sections[nameR.data()];
	logSETUP("Loading config from section %s", nameR.data());

	get_value(section, "n", n);
	get_value(section, "dt", dt);
}
*/
//old routins
/*
template <MODE M, class SEQ, class WF, class PROP, class TASKS, OPTIMS opt>
struct ROUTINE
{
	static constexpr MODE mode = M;
	static constexpr OPTIMS optims = opt;
	using Propagator_t = PROP;
	using Coupling = typename Propagator_t::Coupling;
	// using Source_t = SOURCE;
	using Passes_t = SEQ;
	using Tasks_t = TASKS;
	using WF_t = WF;
	// using COMPS = typename mfilter_all<is_comp, Args...>::type;
	// using DUMPS = typename mfilter_all<is_dump, Args...>::type;
	// using DUMPS = typename filter<is_dump, tuple, Args...>::type;
	string_view name;
	DATA_FORMAT data_format;
	DUMP_FORMAT dump_format;
	ADV_CONFIG advanced;
	constexpr explicit ROUTINE(string_view name,
							   DATA_FORMAT data_format,
							   DUMP_FORMAT dump_format,
							   ADV_CONFIG advanced
	) : name(name),
		data_format(data_format),
		dump_format(dump_format),
		advanced(advanced) {}

	// inline bool load(WF_t& wf)
	// {
	// 	if constexpr (is_same_v<Source_t, no_source>)
	// 		return false;
	// 	else
	// 	{
	// 		wf.load(Source_t);
	// 		return true;
	// 	}
	// }
	inline bool stepExitCondition();
	inline bool stateExitCondition();


	constexpr auto Mode() const { return mode; }
	constexpr auto Name() const { return name; }
	constexpr auto Optims() const { return optims; }
	// constexpr auto startsInP() const { return advanced.propagator.startsInP(); }
	constexpr auto discardPhase() const { return advanced.discard_eigenstate_phase; }
	constexpr auto PropagatorSplits() { return Propagator_t::splits; }

	bool canRun() const { return ((mode & MODES) || MODES == 0); }
	void run() const
	{
		if ((mode & MODES) || MODES == 0)
		{
			logInfo("Routine %s will be run", name.data());
		}
		else
		{
			logWarning("Routine %s will not be run", name.data());
		}
	}
};

template <typename T, typename PROP, typename SEQ, typename...Args>
struct ROUTINE
{
	using Propagator_t = PROP;
	using Regions_t = typename Propagator_t::Regions;
	// using COMPS = typename mfilter_all<is_comp, Args...>::type;
	MODE mode;
	string_view name;
	SEQ passes;
	OPTIMS optims;
	DATA_FORMAT data_format;
	DUMP_FORMAT dump_format;
	PROP propagator;
	T advanced;
	tuple<Args...> tasks;
	constexpr explicit ROUTINE(MODE mode,
							   string_view name,
							   SEQ passes,
							   OPTIMS optims,
							   PROP propagator,
							   DATA_FORMAT data_format,
							   DUMP_FORMAT dump_format,
							   T advanced,
							   Args... args
	) : mode(mode),
		name(name),
		passes(passes),
		optims(optims),
		data_format(data_format),
		dump_format(dump_format),
		propagator(propagator),
		advanced(advanced),
		tasks{ args... } { }


	constexpr bool imaginaryTimeQ() const { return mode == IM; }
	inline bool stepExitCondition();
	inline bool stateExitCondition();
	constexpr auto Mode() const { return mode; }
	constexpr auto Name() const { return name; }
	constexpr auto Optims() const { return optims; }
	constexpr auto propagatorName() const { return propagator.name; }
	// constexpr auto startsInP() const { return advanced.propagator.startsInP(); }
	constexpr auto disablesCAP() const { return advanced.CAP_opt.disable; }
	constexpr auto CAPlength() const { return advanced.CAP_opt.length_au; }
	constexpr auto discardPhase() const { return advanced.discard_eigenstate_phase; }
	constexpr auto PropagatorSplits() { return advanced.propagator.splits; }
	constexpr auto SplitPropagator() { return propagator; }

	bool canRun() const { return ((mode & MODES) || MODES == 0); }
	void run() const
	{
		if ((mode & MODES) || MODES == 0)
		{
			logInfo("Routine %s will be run", name.data());
		}
		else
		{
			logWarning("Routine %s will not be run", name.data());
		}
	}
};

	, AFTER < DUMP_PSI < REP::X, DIMflag>>
{
	join_v <STD_VKV_PROP<H>::name>,
		DATA_FORMAT{ true },
		DUMP_FORMAT{ true, true, true },
		ADV_CONFIG
	{
		.load_from_source = false,
		.discard_eigenstate_phase = true,
		.fs_postpropagation = 0.0
	}
};

template <size_t PASSES, typename WF, typename H, OPTIMS opt>
constexpr auto STANDARD_IM()
{
	return ROUTINE<IM, up_to_t<PASSES>, WF, STD_VKV_PROP<H>, STD_IM<PASSES, H>, opt>
	{
		join_v <STD_VKV_PROP<H>::name>,
			DATA_FORMAT{ true },
			DUMP_FORMAT{ true, true, true },
			ADV_CONFIG
		{
			.load_from_source = false,
			.discard_eigenstate_phase = true,
			.fs_postpropagation = 0.0
		}
	};
}

template <size_t PASSES, typename WF, typename H, OPTIMS opt>
constexpr auto STANDARD_RE()
{
	return ROUTINE<RE, up_to_t<PASSES>, WF, STD_VKV_PROP<H>, STD_RE<PASSES, H>, opt>
	{
		join_v <STD_VKV_PROP<H>::name>,
			DATA_FORMAT{ .binary = true },
			DUMP_FORMAT{ .binary = true, .complex = false, .unnormalized = false },
			ADV_CONFIG
		{
			.load_from_source = true,
			.discard_eigenstate_phase = true,
			.fs_postpropagation = 0.0
		}
	};
}



template <size_t PASSES, typename PROP>
constexpr auto STANDARD_IM(PROP prop)
{
	using Coupling = typename PROP::Coupling;
	return ROUTINE{
		IM,
		join_v <PROP::name>,
		from_to<0, PASSES>,
		NO_OPTIMIZATIONS,
		prop,
		DATA_FORMAT{ true },
		DUMP_FORMAT{ true, true, true },
		ADV_CONFIG
		{
			.load_from_source = false,
			.discard_eigenstate_phase = true,
			.fs_postpropagation = 0.0
		},
		DIRECT_VALUES<StepOperator>(),
		CO_ORTHONORMALIZE_I<from_to_t<0, PASSES>>(),
		IMMEDIATE_AVG<typename PROP::Kinetic>(),
		IMMEDIATE_AVG<typename PROP::Potential>(),
		DIRECT_VALUES<EnergySumOperator, EnergyDiffOperator>(),
		DUMP_PSI<AFTER>{AFTER{}, DIMflag, REP::X, 1000000.0}
	};
}

template <size_t PASSES, typename PROP>
constexpr auto STANDARD_RE(PROP prop)
{
	using Coupling = typename PROP::Coupling;
	return ROUTINE
	{
		RE,
		join_v <PROP::name>,
		from_to<0, PASSES>,
		NO_OPTIMIZATIONS,
		prop,
		DATA_FORMAT{.binary = true },
		DUMP_FORMAT{.binary = true, .complex = false, .unnormalized = false },
		ADV_CONFIG
		{
			.load_from_source = true,
			.discard_eigenstate_phase = true,
			.fs_postpropagation = 0.0
		}
		, DIRECT_VALUES<StepOperator, TimeOperator>()
		, Coupling()
		// , CO_AVG < typename PROP::Kinetic>()
		// , CO_AVG < typename PROP::Potential>()
		, CO_AVG < IdentityOperator<>>()
		// , CO_DIR_AVG <PotentialDerivativeX<typename PROP::Potential>>()
		// , CO_PROJ<IdentityOperator<>, seq<0>>()
		// , CO_FLUXES<PROP
		// ,BOX<10>
		// ,N2S, S2D, N2D
		// , S2CAP
		// , D2CAP
		// , T2CAP
		// >()
		, AUXILLARY_VALUES<ETAOperator>()
		, DUMP_PSI(AFTER{}, D3, REP::X)
		, DUMP_PSI(BEFORE{}, D3, REP::X)
		, DUMP_PSI(DURING{.step = 1 }, D3, REP::X)
		, DUMP_PSI(DURING{.interval = 50 }, D3, REP::X)
	};
}
*/