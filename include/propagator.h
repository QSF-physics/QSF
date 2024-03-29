#include <chrono>
#include "splitting.h"

namespace QSF
{
#pragma region AvailableOperations
	struct Step
	{
		static constexpr REP rep	= REP::NONE;
		static constexpr bool late= false;
	};
	struct Time
	{
		static constexpr REP rep	= REP::NONE;
		static constexpr bool late= false;
	};

	struct ETA
	{
		static constexpr REP rep										= REP::NONE;
		static constexpr bool late									= false;
		static constexpr std::string_view name			= "ETA [h:m:s]";
		static constexpr std::string_view formatting= "%8d:%02d:%02d";
	};

	using ENERGY_TOTAL		 = SUM<AVG<PotentialEnergy>, AVG<KineticEnergy>>;
	using ENERGY_DIFFERENCE= CHANGE<SUM<AVG<PotentialEnergy>, AVG<KineticEnergy>>>;

#pragma endregion AvailableOperations
	struct PropagatorBase
	{
		double dt;
		ind max_steps;
		double state_accuracy;

		ind step{0};
		double timer{0.0};
		double timer_copy{0.0};

		inline void incrementBy(double fraction= 1.0) { timer+= dt * fraction; }
		inline void time_backup() { timer_copy= timer; }
		void time_restore() { timer= timer_copy; }
		void time_reset() { timer= 0; }
		void reset()
		{
			step = 0;
			timer= 0;
		}
	};

	template<MODE M, class SpType, class HamWF> struct SplitPropagator: Config, PropagatorBase
	{
		using clock	 = std::chrono::high_resolution_clock;
		using TimeVar= clock::time_point;

		TimeVar pass_start;

		using SplitType		 = SpType;
		using ChainExpander= typename SplitType::ChainExpander;

		template<uind chain> using Chain= typename SplitType::template Chain<chain>;

		template<uind chain> using splits= typename Chain<chain>::splits;

		static constexpr uind ChainCount			= ChainExpander::size;
		static constexpr REP firstREP					= SplitType::firstREP;
		static constexpr REP invREP						= REP::BOTH ^ firstREP;
		static constexpr MODE mode						= M;
		static constexpr std::string_view name= modeName(M);	 // == MODE::IM ? "IM" : "RE";

		// static constexpr std::string_view name = SplitType::name;
		Section settings;
		HamWF wf;
		double tot_energy;
		double dif_energy;
		cxd* psi_copy= nullptr;		// backup ψ(x,t) on m_l grid
		cxd* psi_acc = nullptr;		// accumulated ψ(x,t) on m_l grid

#pragma region Auto
		void autosetTimestep()
		{
			if(dt == 0.0)
			{
				// dt = 0.5 * pow(0.5, DIM - 1) * dx * dx;
				logWarning("Timestep dt was 0, using automatic value dt=%g", dt);
			}
		}
		// TODO
		//  IO::path save(IO::path path = "", DUMP_FORMAT df = { .dim = DIM, .rep = firstREP })
		//  {

		// }
		// void load(IO::path input_path)
		// {
		// 	if (!input_path.has_extension())
		// 		input_path += IO::psi_ext;
		// 	wf._load(input_path);

		// }
#pragma endregion Auto

#pragma region Initialization
		/// @brief Constructs a propagator
		/// @tparam M Type of the evolution mode mode (MODE::IM or MODE::RE)
		/// @tparam SpType Type of the split operator
		/// @tparam HamWF Type of the Hamiltonian implicitly defined wavefunction
		/// @param pb Propagator base instance or its initializer-list
		/// @param wf Wavefunction (and hence grid) on which to operate on
		/// @param source_section Section of the ini file from which to read config
		// (defaults to mode name)
		explicit SplitPropagator(
			PropagatorBase pb, HamWF&& wf, std::string source_section= std::string(name))
			: Config(source_section),
				PropagatorBase(pb),
				wf(wf)
		{
			logInfo("SplitPropagator init");
			if constexpr(M == MODE::RE)
			{
				max_steps= wf._coupling.maxPulseDuration() / dt + 1;
				logSETUP(
					"maxPulseDuration: %g, dt: %g => max_steps: %td",
					wf._coupling.maxPulseDuration(),
					dt,
					max_steps);
				state_accuracy= 0;
			}
			init();
		}
		// SplitPropagator(HamWF&& wf, std::string source_section = std::string(name))
		// 	: Config(source_section), wf(wf)
		// {
		// 	logInfo("SplitPropagator init");
		// 	max_steps = wf._coupling.maxPulseDuration() / dt + 1;
		// 	logSETUP("maxPulseDuration: %g, dt: %g => max_steps: %td", wf._coupling.maxPulseDuration(),
		// dt, max_steps); 	state_accuracy = 0; 	init();
		// }
		SplitPropagator(
			std::string source_section= std::string(name), std::string config_file= "project.ini")
			: Config(source_section, config_file),
				settings(ini.sections[source_section]),
				wf(settings)
		{
			inipp::get_value(settings, "dt", dt);
			if constexpr(M == MODE::IM)
			{
				inipp::get_value(settings, "max_steps", max_steps);
				inipp::get_value(settings, "state_accuracy", state_accuracy);
				logInfo("max_steps %td state_accuracy %g", max_steps, state_accuracy);
			}
			else
			{
				max_steps= wf._coupling.maxPulseDuration() / dt + 1;
				logSETUP(
					"maxPulseDuration: %g, dt: %g => max_steps: %td",
					wf._coupling.maxPulseDuration(),
					dt,
					max_steps);
			}
			init();
		}
		void init()
		{
			if constexpr(ChainCount > 1)
			{
				logInfo("Using Operator Split Groups (Multiproduct splitting)");
				psi_copy= new cxd[wf.localSize()];
				psi_acc = new cxd[wf.localSize()];
			}

			file_log= IO::fopen_with_check((source_section + IO::log_ext), "w");
			logInfo("dt %g", dt);
		}
		~SplitPropagator()
		{
			if constexpr(ChainCount > 1)
			{
				logInfo("Using Operator Split Groups (Multiproduct splitting)");
				delete[] psi_copy;
				delete[] psi_acc;
				psi_copy= nullptr;
				psi_acc = nullptr;
			}
			logInfo("SplitPropagation done!\n");
		}
#pragma endregion Initialization

#pragma region Computations
		template<class Op> inline double getValue()
		{
			// LOG_INLINE(formatting.data(), val / 3600, val % 3600 / 60, val % 60);
			if constexpr(std::is_same_v<Time, Op>) return timer;
			if constexpr(std::is_same_v<Step, Op>) return step;
			if constexpr(std::is_same_v<ETA, Op>)
			{
				using std::chrono::duration;
				using std::chrono::duration_cast;
				using std::chrono::milliseconds;
				using std::chrono::seconds;
				auto sec= duration_cast<seconds>(clock::now() - pass_start).count();
				// TODO: use atstep
				return (1.0 * max_steps / step - 1.0) * sec;
			}
			else if constexpr(!std::is_same_v<Step, Op> && !std::is_same_v<Time, Op>)
				return wf.template getValue<Op>();
		}
		template<REP R, class BO, class... Op> inline void compute(BO& bo, VALUE<Op...>&&)
		{
			using T= VALUE<Op...>;
			bo.template store<T>(getValue<Op>()...);
		}

		template<REP R, class BO, class... Op> inline void compute(BO& bo, SUM<Op...>&&)
		{
			using T= SUM<Op...>;
			bo.template store<T>((bo.template getLastValue<Op>() + ...));
			if constexpr(std::is_same_v<T, ENERGY_TOTAL>) tot_energy= bo.template getLastValue<T>();
		}

		template<REP R, class BO, class Op> inline void compute(BO& bo, CHANGE<Op>&&)
		{
			static double last= 0;
			using T						= CHANGE<Op>;
			double curr				= bo.template getLastValue<Op>();
			MPI::reduceImmediataly(&curr);
			bo.template store<T>(curr - last);
			if constexpr(std::is_same_v<T, ENERGY_DIFFERENCE>) dif_energy= curr - last;
			last= curr;
		}

		template<REP R, class BO, class... Op> inline void compute(BO& bo, OPERATION<Op...>&& c)
		{
			((wf.template operation<M, R, Op>()), ...);
		}

		template<REP R, class BO, class... Op> inline void compute(BO& bo, AVG<Op...>&& c)
		{
			using T= AVG<Op...>;
			bo.template store<T>((wf.template average<R, Op>())...);
		}

		template<REP R, class BO, uind dir, class Op>
		inline void compute(BO& bo, AVG<DERIVATIVE<dir, Op>>&& c)
		{
			using T= AVG<DERIVATIVE<dir, Op>>;
			bo.template store<T>((wf.template average_der<R, dir, Op>()));
		}

		template<uind... Is> inline void velocityGaugeCorrection(borVec<HamWF::DIM>* curr, seq<Is...>)
		{
			ind counter= 0;
			while(counter < wf.m_l)
			{
				((curr[counter][Is]+= wf._coupling[Axis<Is>] * norm(wf[counter])), ...);
				counter++;
			}
		}

		template<REP R, class BO, class... Op> inline void compute(BO& bo, FLUX<Op...>&& c)
		{
			using T											= FLUX<Op...>;
			borVec<HamWF::DIM>* currents= wf.template current_map();
			if constexpr(HamWF::couplesInRep == REP::P) velocityGaugeCorrection(currents, HamWF::indices);
			bo.template store<T>((wf.template flux<R, Op>())...);
		}

		// If no match here is found pass to the wavefunction with buffer
		template<REP R, class BO, class COMP> inline void compute(BO& bo, COMP&& c)
		{
			wf.template compute<M, R>(bo, std::forward<COMP>(c));
		}

		inline void ditch() {}

		template<WHEN when, bool B, class... COMP>
		inline void computeEach(BufferedOutputs<B, COMP...>& bo)
		{
			if(!HamWF::MPIGridComm::calcsEnabled()) return;
			// Run only if required REP is firstREP or no preference
			((!COMP::late && (COMP::rep == REP::NONE || bool(COMP::rep & firstREP))
					? compute<firstREP>(bo, COMP{})
					: ditch()),
			 ...);

			// Check if work in inverse fourier space is needed
			if constexpr((false || ... || (bool(COMP::rep & invREP))))
			{
				wf.template fourier<invREP>();
				((bool(COMP::rep & invREP) ? compute<invREP>(bo, COMP{}) : ditch()), ...);
				wf.template fourier<firstREP>();
			}

			// Check if any "late" operations are required
			if constexpr((false || ... || (COMP::late)))
			{
				((COMP::late ? compute<invREP>(bo, COMP{}) : ditch()), ...);
			}

			bo.template logOrPass<M, when>(step);
		}
#pragma endregion Computations

#pragma region Evolution
		inline bool shouldContinue()
		{
			if constexpr(mode == MODE::RE) return (step <= max_steps);
			else
				return fabs(dif_energy) > state_accuracy && step != max_steps;
		}

		template<uind chain, uind... SI> inline void chainEvolve(seq<SI...>)
		{
			(
				[&]
				{
					constexpr REP rep= Chain<chain>::template rep<SI>;
					if(SI > 0) wf.template fourier<rep>();

					if constexpr(REP::BOTH == HamWF::couplesInRep || rep == HamWF::couplesInRep)
						wf.template precalc<rep, OPTIMS::NONE>(timer);
					//  Timings::measure::start(op.name);
					wf.template evolve<M, rep>(dt * Chain<chain>::mults[SI]);
					// Timings::measure::stop(op.name);
					// if (step == 4)
					// printf("SplitGroup %td Evolving in REP %td with delta=%g dt=%g pID:%d\n", chain,
					// ind(rep), Chain<chain>::mults[SI], dt, MPI::pID);
					if constexpr(REP::BOTH == HamWF::couplesInRep) incrementBy(Chain<chain>::mults[SI] * 0.5);
					else if constexpr(rep == HamWF::couplesInRep)
						incrementBy(Chain<chain>::mults[SI] * 1.0);
				}(),
				...);
		}

		template<uind... chain> inline void makeStep(seq<chain...>)
		{
			chainBackup();
			((chainRestore<chain>()
				//   , Timings::measure::start("EVOLUTION")
				,
				chainEvolve<chain>(splits<chain>{})
				//   , Timings::measure::stop("EVOLUTION")
				,
				chainComplete<chain>()),
			 ...);
			step++;
		}

		template<uind chain, uind... SI> inline void chainFakeEvolve(seq<SI...>)
		{
			(
				[&]
				{
					constexpr REP rep= Chain<chain>::template rep<SI>;
					if constexpr(REP::BOTH == HamWF::couplesInRep || rep == HamWF::couplesInRep)
						wf.template precalc<rep, OPTIMS::NONE>(timer);
					if constexpr(REP::BOTH == HamWF::couplesInRep) incrementBy(Chain<chain>::mults[SI] * 0.5);
					else if constexpr(rep == HamWF::couplesInRep)
						incrementBy(Chain<chain>::mults[SI] * 1.0);
				}(),
				...);
		}
		template<uind... chain> inline void makeFakeStep(seq<chain...>)
		{
			((chainFakeEvolve<chain>(splits<chain>{})), ...);
			step++;
		}
#pragma endregion Evolution

#pragma region SplitChains

		template<ind chain> void chainRestore()
		{
			if constexpr(chain > 1)
			{
				if(step == 4) { logInfo("group %td restore", chain); }
				time_restore();
				for(ind i= 0; i < wf.localSize(); i++) wf[i]= psi_copy[i];
				// HamWF::Coupling::restore();
			}
		}
		template<ind chain> void chainComplete()
		{
			if constexpr(ChainCount > 1)
			{
				if(step == 4)
				{
					logInfo("group %td complete, applying coeff=%g", chain, Chain<chain>::value);
				}
				if constexpr(chain == ChainCount - 1)
					for(ind i= 0; i < wf.localSize(); i++) wf[i]= Chain<chain>::value * wf[i] + psi_acc[i];

				else
					for(ind i= 0; i < wf.localSize(); i++) psi_acc[i]+= Chain<chain>::value * wf[i];
			}
		}
		void chainBackup()
		{
			if constexpr(ChainCount > 1)
			{
				if(step == 4) { logInfo("group backup because ChainCount=%td", ChainCount); }
				for(ind i= 0; i < wf.localSize(); i++)
				{
					psi_copy[i]= wf[i];
					psi_acc[i] = 0.0;
				}
				time_backup();
				// HamWF::Coupling::backup();
			}
		}
#pragma endregion splitChains

#pragma region Runner

		/// @brief
		/// @tparam OUTS (auto-deduced)
		/// @tparam Worker (auto-deduced)
		/// @param outputs Instance of BufferedOutputs (operations to be performed each step)
		/// @param worker Function of args (const WHEN when, const ind step, const uind pass, auto& wf)
		/// @param restore Whether to restore the calculations from the last snapshot
		template<class OUTS, class Worker> void run(OUTS&& outputs, Worker&& worker, bool restore= 0)
		{
			pass_start= clock::now();

			ind atstep= 0;
			if(restore)
			{
				FILE* f= IO::fopen_with_check("backup.info", "rb");
				fread(&atstep, sizeof(ind), 1, f);
				fclose(f);
				printf("pID %d, Read step: %td\n", MPI::pID, atstep);
				outputs.template init<M>(source_section, atstep);
				// TODO: make this for gauges interacting with potentials (not fields)
				while(step < atstep) makeFakeStep(ChainExpander{});
				logWarning("Fast forwarded computation to step %td", step);
				wf.restore();
			}
			else
			{
				outputs.template init<M>(source_section, 0);
				worker(WHEN::AT_START, step, 0, wf);
				computeEach<WHEN::AT_START>(outputs);
			}
			while(shouldContinue())
			{
				makeStep(ChainExpander{});
				computeEach<WHEN::DURING>(outputs);
				wf.postCompute();
				// logWarning("%d", HamWF::MPIGridComm::many);
				worker(WHEN::DURING, step, 0, wf);
			}
			makeStep(ChainExpander{});
			computeEach<WHEN::AT_END>(outputs);
			wf.postCompute();
			worker(WHEN::AT_END, step, 0, wf);

			if(!MPI::pID)
			{
				std::filesystem::remove("backup.info");
				logInfo("backup.info file removed");
			}
			wf.remove_backups();
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

		template<class OUTS> void run(OUTS&& outputs, bool restore= 0)
		{
			run<OUTS>(
				std::forward<OUTS>(outputs),
				[](WHEN when, ind step, uind pass, const auto& wf) {},
				restore);
		}

		template<class OUTS, class Worker> void run(Worker&& worker, bool restore= 0)
		{
			run(std::move(OUTS{settings}), std::forward<Worker>(worker), restore);
		}

		template<class OUTS> void run(bool restore= 0)
		{
			run<OUTS>([](WHEN when, ind step, uind pass, const auto& wf) {}, restore);
		}

#pragma endregion Runner
	};

	template<MODE M, class SpType, class HamWF> auto SplitPropagate(HamWF wf)
	{
		return SplitPropagator<M, SpType, HamWF>(wf);
	}
};	 // namespace QSF