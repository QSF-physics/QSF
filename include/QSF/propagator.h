#include "splitting.h"

struct PropagatorBase
{
	double dt;
	ind step{ 0 };
	double timer{ 0.0 };

	PropagatorBase(Section& settings)
	{
		inipp::get_value(settings, "dt", dt);
		logInfo("dt %g", dt);
		//TODO: if 0 autoset?
	}
	void reset()
	{
		step = 0; timer = 0;
	}
	void incrementBy(double fraction = 1.0) { timer += dt * fraction; }
};

template <MODE M>
struct SplitPropagatorBase;

template <>
struct SplitPropagatorBase<MODE::IM> : PropagatorBase
{
	ind max_steps;
	double state_accuracy;
	double dE = 10.0;
	double energy;
	bool stillEvolving()
	{
		return fabs(dE) > state_accuracy && step != max_steps;
		// !((fabs(dE) < state_accuracy || dE / energy == 0.) || step == max_steps);
	}
	SplitPropagatorBase(Section& settings) : PropagatorBase(settings)
	{
		inipp::get_value(settings, "max_steps", max_steps);
		inipp::get_value(settings, "state_accuracy", state_accuracy);
		logInfo("max_steps %td", max_steps);
	}
};

template <>
struct SplitPropagatorBase<MODE::RE> : PropagatorBase
{
	ind max_steps;
	bool stillEvolving()
	{
		return (step < max_steps);
	}
	SplitPropagatorBase(Section& settings) : PropagatorBase(settings)
	{
		// inipp::get_value(settings, "max_steps", max_steps);
	}
};


struct Step {
	static constexpr REP rep = REP::BOTH;
};
struct Time {
	static constexpr REP rep = REP::BOTH;
};

template <MODE M, class HamWF, class SpBase, template <class, size_t> class SpType, size_t Order>
struct SplitPropagator : SplitPropagatorBase<M>
{
	using SplitPropagatorBase<M>::step;
	using SplitPropagatorBase<M>::timer;
	using SplitPropagatorBase<M>::dt;
	using SplitPropagatorBase<M>::incrementBy;
	using SplitPropagatorBase<M>::reset;
	using SplitBase = SpBase;
	using SplitType = SpType<SplitBase, Order>;
	using ChainExpander = typename SplitType::ChainExpander;

	template <size_t splitGroup>
	using Chain = typename SplitType::template Chain<splitGroup>;

	template <size_t splitGroup>
	using reps = typename Chain<splitGroup>::reps;
	template <size_t splitGroup>
	using splits = typename Chain<splitGroup>::splits;

	static constexpr size_t ChainCount = ChainExpander::size;
	static constexpr REP startsIn = SplitBase::firstREP;
	// static constexpr REP couplesInRep = C::couplesInRep;
	static constexpr MODE mode = M;
	static constexpr std::string_view name = (M == MODE::IM) ? "IM" : "RE";
	// static constexpr std::string_view name = SplitType::name;
	HamWF wf;
	SplitPropagator(Section& settings) : SplitPropagatorBase<M>(settings), wf(settings)
	{

	}


	double timer_copy;
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
	template <REP R, class Op>
	inline double calc()
	{
		if constexpr (std::is_same_v<Op, Time>)
			return timer;
		else if constexpr (std::is_same_v<Op, Step>)
			return double(step);
		else return wf.template avg<AXIS::XYZ, Op>();
	}

	template <size_t... splitGroup>
	void makeStep(seq<splitGroup...> t)
	{
		step++;
		// (printf("split %td %td\n", splitGroup, ChainCount), ...);
		groupBackup();
		((
			groupRestore<splitGroup>()
		//   , Timings::measure::start("EVOLUTION")
			, evolve<splitGroup>(reps<splitGroup>{}, splits<splitGroup>{})
		  //   , Timings::measure::stop("EVOLUTION")
			, groupComplete<splitGroup>()
			), ...);
		  // wf.post_evolve();
	}

	static bool calcsEnabled()
	{
		return true;
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