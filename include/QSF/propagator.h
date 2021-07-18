#include "splitting.h"

struct Step {
	static constexpr REP rep = REP::BOTH;
};
struct Time {
	static constexpr REP rep = REP::BOTH;
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

	template <uind chain>
	using Chain = typename SplitType::template Chain<chain>;
	template <uind chain>
	using reps = typename Chain<chain>::reps;
	template <uind chain>
	using splits = typename Chain<chain>::splits;

	static constexpr uind ChainCount = ChainExpander::size;
	static constexpr REP firstREP = SplitType::firstREP;
	static constexpr std::string_view name = M == MODE::IM ? "IM" : "RE";
	static constexpr MODE mode = M;
	// static constexpr REP couplesInRep = C::couplesInRep;

	Section settings;
	// static constexpr std::string_view name = SplitType::name;
	HamWF wf;
	double energy{ 0 };
	double dE{ 10.0 };
	// SplitPropagator(Section& settings) : PropagatorBase(settings), TimeMode<M>(settings), wf(settings), outputs(settings, M) {}
	void autosetTimestep()
	{
		if (dt == 0.0)
		{
			// dt = 0.5 * pow(0.5, DIM - 1) * dx * dx;
			logWarning("Timestep dt was 0, using automatic value dt=%g", dt);
		}
	}

	SplitPropagator(PropagatorBase pb, HamWF wf) : PropagatorBase(pb), wf(wf) {
		max_steps = 100;//FIXME 
		state_accuracy = 0;
	}

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
		if constexpr (ChainCount > 1) wf.initHelpers();

		file_log = openLog(name);
		logInfo("dt %g", dt);
	}

	~SplitPropagator()
	{
		logInfo("SplitPropagation done!\n");
	}



	// template <typename Op> inline void getOperator(Op);
	inline double getOperator(Time) { return timer; }
	inline double getOperator(Step) { return step; }

	//If no match here is found pass to the wavefunction
	template < REP R, class BO, class COMP, uind...Is>
	inline void compute(BO& bo, COMP&& c, seq<Is...>&& s)
	{
		wf.template compute<R>(bo, std::move(c), std::move(s));
	}
	template <REP R, class BO, class... Op, uind...Is>
	inline void compute(BO& bo, PROPAGATOR_VALUE<Op...>&&, seq<Is...>&&)
	{
		using T = PROPAGATOR_VALUE<Op...>;
		// logInfo("returning pos %td %td", Is...);
		((bo.template storeInBuffer < Is, T>(getOperator(Op{}))), ...);
	}

	inline void ditch() {}

	template <bool B, class... COMP>
	inline void computeEach(BufferedOutputs<B, COMP...>& bo)
	{
		using BO = BufferedOutputs<B, COMP...>;

		((COMP::template canRun<firstREP, EARLY>
		  ? compute<firstREP>(bo, COMP{}, BO::template pos_v<COMP>())
		  : ditch()),
		 ...);

		if constexpr (BO::template needsFFT<firstREP>())
		{
			wf.template fourier<invREP(firstREP)>();
			((COMP::template canRun<invREP(firstREP), LATE>
			  ? compute<invREP(firstREP)>(bo, COMP{}, BO::template pos_v<COMP>())
			  : ditch()),
			 ...);
			wf.template fourier<firstREP>();
		}
		else if (step < 2)
		{
			logWarning("Computations do not require FFT, if you need to output wavefunction in the opposite REP make sure to export the WF explicitly in the desired REP.");
		}

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

	template <uind chain, uind ... repI, uind ... SI>
	inline void evolve(seq<repI...>, seq<SI...>)
	{
		([&]
		 {
			 constexpr REP rep = REP(repI);
			 wf.template fourier<rep>();
			//  wf.template precalc<rep, NO_OPTIMIZATIONS>();
			//  Timings::measure::start(op.name);
			 wf.template evolve<M, rep>(dt * Chain<chain>::mults[SI]);
			//  Timings::measure::stop(op.name);
			 if (step == 4)
				 logInfo("SplitGroup %td Evolving in REP %td with delta=%g", chain, repI, Chain<chain>::mults[SI]);
			 if constexpr (REP::BOTH == HamWF::couplesInRep)
				 incrementBy(Chain<chain>::mults[SI] * 0.5);
			 else if constexpr (rep == HamWF::couplesInRep)
				 incrementBy(Chain<chain>::mults[SI] * 1.0);
		 }(), ...);
	}

	template <class OUTS, class Worker>
	void run(OUTS&& outputs, Worker&& worker, uind pass = 0)
	{
		outputs.initLogger(M, pass, name);
		computeEach(outputs);
		outputs.template logOrPass<WHEN::AT_START>(step);
		worker(step, pass, wf);
		while (stillEvolving())
		{
			makeStep(ChainExpander{});
			wf.post_step();
			computeEach(outputs);
			outputs.template logOrPass<WHEN::DURING>(step);
			worker(step, pass, wf);
		   // outputs.template logOrPass<DURING<>>(step);
		}
		makeStep(ChainExpander{});
		wf.post_step();
		computeEach(outputs);
		outputs.template logOrPass<WHEN::AT_END>(step);
		worker(step, pass, wf);
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
	template <class OUTS>
	void run(OUTS&& outputs, uind pass = 0)
	{
		run<OUTS>(std::forward<OUTS>(outputs), [](ind step, uind pass, const auto& wf) {}, pass);
	}

	template <class OUTS, class Worker>
	void run(Worker&& worker, uind pass = 0)
	{
		run(std::move(OUTS{ settings }), std::forward<Worker>(worker), pass);
	}

	template <class OUTS>
	void run(uind pass = 0)
	{
		run<OUTS>([](ind step, uind pass, const auto& wf) {}, pass);
	}


	template <uind... chain>
	void makeStep(seq<chain...>)
	{
		groupBackup();
		((
			groupRestore<chain>()
		//   , Timings::measure::start("EVOLUTION")
			, evolve<chain>(reps<chain>{}, splits<chain>{})
		  //   , Timings::measure::stop("EVOLUTION")
			, groupComplete<chain>()
			), ...);
		step++;
	}


	template <ind chain> void groupRestore()
	{
		if constexpr (chain > 1)//changed from 0 14.07
		{
			if (step == 4) { logInfo("group %td restore", chain); }
			time_restore();
			wf.restore();
			// HamWF::Coupling::restore();
		}
	}
	template <ind chain> void groupComplete()
	{
		if constexpr (ChainCount > 1)
		{
			if (step == 4)
			{
				logInfo("group %td complete, applying coeff=%g", chain, Chain<chain>::value);
			}
			if constexpr (chain == ChainCount - 1)
				wf.collect(Chain<chain>::value);
			else wf.accumulate(Chain<chain>::value);

		}
	}
	void groupBackup()
	{
		if constexpr (ChainCount > 1)
		{
			if (step == 4) { logInfo("group backup because ChainCount=%td", ChainCount); }
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
	template <uind Chain>
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