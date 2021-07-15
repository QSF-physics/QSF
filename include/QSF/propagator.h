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
	Config config;
	Section settings;

	double dt;
	ind step{ 0 };
	double timer{ 0.0 };
	double timer_copy;
	PropagatorBase(std::string_view name, int DIMS, int ELEC) :
		config("project.ini", DIMS, ELEC),
		settings(config.ini.sections[name.data()])
	{
		file_log = openLog(name);
		inipp::get_value(settings, "dt", dt);
		logInfo("dt %g", dt);
		//TODO: if 0 autoset?
	}

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

template <MODE M>
struct TimeMode;

template <>
struct TimeMode<MODE::IM>
{
	static constexpr MODE mode = MODE::IM;
	static constexpr std::string_view name = "IM";
	ind max_steps;
	double state_accuracy;
	double dE = 10.0;
	double energy;
	inline bool stillEvolving(ind step)
	{
		return fabs(dE) > state_accuracy && step != max_steps;
		// !((fabs(dE) < state_accuracy || dE / energy == 0.) || step == max_steps);
	}
	TimeMode(Section& settings)
	{
		inipp::get_value(settings, "max_steps", max_steps);
		inipp::get_value(settings, "state_accuracy", state_accuracy);
		logInfo("max_steps %td", max_steps);
	}
};


template <>
struct TimeMode<MODE::RE>
{
	static constexpr MODE mode = MODE::RE;
	static constexpr std::string_view name = "RE";
	ind max_steps;
	bool stillEvolving(ind step)
	{
		return (step <= max_steps);
	}
	TimeMode(Section& settings)
	{
		// inipp::get_value(settings, "max_steps", max_steps);
	}
};

template <MODE M, class SpType, class HamWF>
struct SplitPropagator : PropagatorBase, TimeMode<M>
{
	using SplitType = SpType;
	using ChainExpander = typename SplitType::ChainExpander;

	template <size_t splitGroup>
	using Chain = typename SplitType::template Chain<splitGroup>;

	template <size_t splitGroup>
	using reps = typename Chain<splitGroup>::reps;
	template <size_t splitGroup>
	using splits = typename Chain<splitGroup>::splits;

	static constexpr size_t ChainCount = ChainExpander::size;
	static constexpr REP firstREP = SplitType::firstREP;
	// static constexpr REP couplesInRep = C::couplesInRep;

	// static constexpr std::string_view name = SplitType::name;
	HamWF wf;

	// SplitPropagator(Section& settings) : PropagatorBase(settings), TimeMode<M>(settings), wf(settings), outputs(settings, M) {}


	SplitPropagator() : PropagatorBase(TimeMode<M>::name, 1, 1),
		TimeMode<M>(settings),
		wf(settings) {}

	inline bool stillEvolving()
	{
		return TimeMode<M>::stillEvolving(step);
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

	template <REP R, class Op>
	inline double calc()
	{
		if constexpr (std::is_same_v<Op, Time>)
		{
			return timer;
		}
		else if constexpr (std::is_same_v<Op, Step>)
		{
			return double(step);
		}
		else
		{
			return wf.template avg<AXIS::XYZ, Op>();
		}
	}

	template <class OUTS, class Worker>
	void run(Worker&& worker, uind i = 0)
	{
		OUTS outputs{ settings, M };
		outputs.init(i, TimeMode<M>::name);
		outputs.template compute< M, EARLY, firstREP>(this);
		outputs.template logOrPass<BEFORE<>>(step);
		while (stillEvolving())
		{
			makeStep(ChainExpander{});
			outputs.template logOrPass<DURING<>>(step);
			wf.post_evolve();
			outputs.template compute< M, EARLY, firstREP>(this);
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