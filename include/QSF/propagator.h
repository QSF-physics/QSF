#include "splitting.h"
#include <fstream>
#include "autoconfig.h"

struct Step {
	static constexpr REP rep = REP::BOTH;
};
struct Time {
	static constexpr REP rep = REP::BOTH;
};

struct Config
{
	inipp::Ini<char> ini;
	char config[100];

	explicit Config(std::string_view name, int DIMS, int ELEC)
	{
		logInfo("Parsing %s", name.data());
		// if (MPI::pID)
		std::ifstream is(name.data());
		// ini.clear();
		ini.parse(is);
		// logINI("Raw INI file:");
		// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
		ini.strip_trailing_comments();
		sprintf(config, "%de%dd", DIMS, ELEC);
		// constexpr auto dimelec = STRINGIFY(ELEC) "e" STRINGIFY(DIM) "d";
		ini.default_section(ini.sections[config]);
		ini.default_section(ini.sections["DEFAULT"]);
		ini.interpolate();
		// logINI("Parsed & interpolated project.ini file:");
		// ini.generate(std::cout);
		// if (DEBUG & DEBUG_INI) ini.generate(std::cout);
		// if (DEBUG & DEBUG_INI) MPI_Barrier(MPI_COMM_WORLD);
	}
};

struct PropagatorBase
{
	double dt;
	ind max_steps;
	double state_accuracy;

	ind step{ 0 };
	double timer{ 0.0 };
	double timer_copy{ 0.0 };

	inline void incrementBy(double fraction = 1.0)
	{
		timer += dt * fraction;
	}
	inline void time_backup()
	{
		timer_copy = timer;
	}
	void time_restore()
	{
		timer = timer_copy;
	}
	void time_reset()
	{
		timer = 0;
	}
	void reset()
	{
		step = 0; timer = 0;
	}
};





template <MODE M, class SpType, class HamWF>
struct SplitPropagator : Config, PropagatorBase
{
	using SplitType = SpType;
	using ChainExpander = typename SplitType::ChainExpander;

	template <size_t splitGroup>
	using Chain = typename SplitType::template Chain<splitGroup>;

	template <size_t splitGroup>
	using reps = typename Chain<splitGroup>::reps;
	template <size_t splitGroup>
	using splits = typename Chain<splitGroup>::splits;

	static constexpr MODE mode = M;
	static constexpr std::string_view name = M == MODE::IM ? "IM" : "RE";
	static constexpr size_t ChainCount = ChainExpander::size;
	static constexpr REP firstREP = SplitType::firstREP;
	// static constexpr REP couplesInRep = C::couplesInRep;

	Section settings;
	// static constexpr std::string_view name = SplitType::name;
	HamWF wf;
	double energy{ 0 };
	double dE{ 10.0 };
	// SplitPropagator(Section& settings) : PropagatorBase(settings), TimeMode<M>(settings), wf(settings), outputs(settings, M) {}
	SplitPropagator(PropagatorBase pb, HamWF wf) : PropagatorBase(pb), wf(wf) {}

	SplitPropagator() :
		Config("project.ini", 1, 1),//DIMS, ELEC
		settings(ini.sections[name.data()]),
		wf(settings)
	{
		inipp::get_value(settings, "dt", dt);
		if constexpr (M == MODE::IM)
		{
			inipp::get_value(settings, "max_steps", max_steps);
			inipp::get_value(settings, "state_accuracy", state_accuracy);
			logInfo("max_steps %td", max_steps);
		}
		else
		{
			//TODO: init RE
		}
		file_log = openLog(name);
		logInfo("dt %g", dt);
	}


	inline bool stillEvolving()
	{
		if constexpr (mode == MODE::RE) return (step <= max_steps);
		else return fabs(dE) > state_accuracy && step != max_steps;
	}

	template <REP R>
	void fourier()
	{
		wf.template fourier<R>();
	}

	template <size_t splitGroup, size_t ... repI, size_t ... SI>
	inline void evolve(seq<repI...>, seq<SI...>)
	{
		([&]
		 {
			 constexpr REP rep = REP(repI);
			 wf.template fourier<rep>();
			//  wf.template precalc<rep, NO_OPTIMIZATIONS>();
			//  Timings::measure::start(op.name);
			 wf.template evolve<M, rep>(dt * Chain<splitGroup>::mults[SI]);
			//  Timings::measure::stop(op.name);
			//  printf("Evolving with %g (%g)\n", dt * Chain<splitGroup>::mults[SI], Chain<splitGroup>::mults[SI]);
			 if constexpr (REP::BOTH == HamWF::couplesInRep)
				 incrementBy(Chain<splitGroup>::mults[SI] * 0.5);
			 else if constexpr (rep == HamWF::couplesInRep)
				 incrementBy(Chain<splitGroup>::mults[SI] * 1.0);
		 }(), ...);
	}


	// template <typename Op> inline void getOperator(Op);
	inline double getOperator(Time) { return timer; }
	inline double getOperator(Step) { return step; }

	//If no match here is found pass to the wavefunction
	template < REP R, class BO, class COMP, size_t...Is>
	inline void compute(BO& bo, COMP&& c, seq<Is...>&& s)
	{
		wf.template compute<R>(bo, std::move(c), std::move(s));
	}
	template <REP R, class BO, class... Op, size_t...Is>
	inline void compute(BO& bo, PROPAGATOR_VALUE<Op...>&&, seq<Is...>&&)
	{
		using T = PROPAGATOR_VALUE<Op...>;
		// logInfo("returning pos %td %td", Is...);
		((bo.template storeInBuffer < Is, T>(getOperator(Op{}))), ...);
	}

	inline void ditch() {}

	template <WHEN when, bool B, class... COMP>
	inline void computeEach(BufferedOutputs<B, COMP...>& bo)
	{
		// compute< M, EARLY, firstREP>(this);
		using BO = BufferedOutputs<B, COMP...>;
		((COMP::template canRun<firstREP, EARLY>
		  ? compute<firstREP, BO>(bo, COMP{}, BO::template pos_v<COMP>())
		  : ditch()),
		 ...);

		if constexpr (BO::template needsFFT<firstREP>())
		{
			wf.template fourier<REP::BOTH^ firstREP>();

			wf.template fourier<firstREP>();
		}
		else if (step < 2)
		{
			logWarning("Computations do not require FFT, if you need to output wavefunction in the opposite REP make sure to export the WF explicitly in the desired REP.");
		}
		bo.template logOrPass<when>(step);
	}

	template <class OUTS, class Worker>
	void run(Worker&& worker, uind i = 0)
	{
		OUTS outputs{ settings, M ,i, name };
		computeEach<WHEN::AT_START>(outputs);
		while (stillEvolving())
		{
			makeStep(ChainExpander{});
			wf.post_step();
			computeEach<WHEN::DURING>(outputs);

		   // outputs.template logOrPass<DURING<>>(step);
		}
		// HACK: This makes sure, the steps are always evenly spaced
		// step += (outputs.comp_interval - 1);
		// Evolution::incrementBy(outputs.comp_interval);
		/* if (imaginaryTimeQ)
		{
			Eigen::store<startREP>(state, PASS, energy);
			Eigen::saveEnergyInfo(RT::name, state, PASS, energy, dE);
			energy = 0;
			energy_prev = 0;
			dE = 0;
		}  */
	}

	template <size_t... splitGroup>
	void makeStep(seq<splitGroup...> t)
	{
		groupBackup();
		((
			groupRestore<splitGroup>()
		//   , Timings::measure::start("EVOLUTION")
			, evolve<splitGroup>(reps<splitGroup>{}, splits<splitGroup>{})
		  //   , Timings::measure::stop("EVOLUTION")
			, groupComplete<splitGroup>()
			), ...);
		step++;
	}


	template <ind splitGroup> void groupRestore()
	{
		if constexpr (splitGroup > 1)//changed from 0 14.07
		{
			// if (step == 2)logInfo("group %td restore", splitGroup();
			time_restore();
			wf.restore();
			// HamWF::Coupling::restore();
		}
	}
	template <ind splitGroup> void groupComplete()
	{
		if constexpr (ChainCount > 1)
		{
			if (step == 2) logInfo("group %td complete", splitGroup);
			constexpr auto coeff = Chain<splitGroup>::value;// <splitGroup>().coeff;
			if constexpr (splitGroup == ChainCount - 1)
				wf.collect(coeff);
			else
				wf.accumulate(coeff);
		}
	}
	void groupBackup()
	{
		if constexpr (ChainCount > 1)
		{
			if (step == 2) logInfo("group %td(size) backup", ChainCount);
			wf.backup();
			time_backup();
			// HamWF::Coupling::backup();
		}
	}

};

struct ADV_CONFIG
{
	bool discard_eigenstate_phase = true;
	double fs_postpropagation = 0.0;
	DATA_FORMAT data_format;
	DUMP_FORMAT dump_format;
};

/* Imagi time
struct ImagTime
{
	static constexpr std::string_view name = "IM";
	static constexpr MODE mode = IM;
	double operator()(double val)
	{
		return exp(-val);
	}
	inline bool exitCondition()
	{
		return !((fabs(dE) < state_accuracy || dE / energy == 0.) || step == max_imaginary_steps);
	}
	void config(Section& settings)
	{
		inipp::get_value(settings, "max_imaginary_steps", max_imaginary_steps);
		inipp::get_value(settings, "state_accuracy", state_accuracy);
	}
	inline void after()
	{
		// basePass::after();
		// Eigen::store<startREP>(state, PASS, energy);
		// Eigen::saveEnergyInfo(name, state, PASS, energy, dE);
		// energy = 0;
		// energy_prev = 0;
		// dE = 0;
	}
};
struct RealTime
{
	static constexpr std::string_view name = "RE";
	static constexpr MODE mode = RE;
	cxd operator()(double val)
	{
		return cos(-val) + I * sin(-val);
	}
	inline bool exitCondition()
	{
		return (step < ntsteps);
	}
};
*/

/* old get chain
SplitType splitSum;
	template <size_t Chain>
	constexpr size_t splitCount() const { return splitSum.sizes[Chain]; };

	constexpr SplitPropagator() {}
	template <typename ... deltaMult>
	constexpr SplitPropagator(SplitType splitSum, double CAPlength, deltaMult...mult) :
		splitSum(splitSum) {}
			splits{ SplitOperators<T_Op, V_Op, C_Op>{mult, CAPlength} ... } {}


	template <size_t chain> constexpr auto getChain() const { return get<chain>(splitSum.chains); }

	template <size_t chain, size_t I> constexpr auto getOperator() const
	{
		return  get<I>(getChain<chain>().splits);
	}

	*/