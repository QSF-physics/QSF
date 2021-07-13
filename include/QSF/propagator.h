
#include "splitting.h"

struct PropagatorBase
{
	double timer_copy;
	double timer;
	double dt;
	ind step;
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
	void incrementBy(double fraction) { timer += dt * fraction; }
};

template <MODE M, class HamWF, class SpBase, template <class, size_t> class SpType, size_t Order>
struct SplitPropagator : PropagatorBase//TODO: rename to SplitPropagator
{
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
	static constexpr std::string_view name = M == MODE::IM ? "IM" : "RE";
	// static constexpr std::string_view name = SplitType::name;
	HamWF wf;
	SplitPropagator() {}

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

			 if constexpr (REP::BOTH == HamWF::couplesInRep) time_incrementBy(Chain<splitGroup>::mults[SI] * 0.5);
			 else if constexpr (rep == HamWF::couplesInRep) time_incrementBy(Chain<splitGroup>::mults[SI] * 1.0);
		 }(), ...);
	}

	template <ind splitGroup> void groupRestore()
	{
		if constexpr (splitGroup > 0)
		{
			// if (step == 2)logInfo("group %td restore", splitGroup();
			time_restore();
			wf.restore();
			HamWF::Coupling::restore();
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
			HamWF::Coupling::backup();
		}
	}


	template <size_t... splitGroup>
	void makeStep(seq<splitGroup...> t)
	{
		groupBackup(wf);
		((groupRestore<splitGroup>(wf)
		//   , Timings::measure::start("EVOLUTION")
		  , evolve<splitGroup>(reps<splitGroup>{}, splits<splitGroup>{}, wf)
		//   , Timings::measure::stop("EVOLUTION")
		  , groupComplete<splitGroup>(wf)
		  ), ...);
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