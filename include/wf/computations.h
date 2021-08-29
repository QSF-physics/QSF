struct _COMPUTATION {};
struct _ALLOCS_EIGENSTATES {};
// template <typename SEQ> struct ALLOCS_EIGENSTATES : _ALLOCS_EIGENSTATES
// {
// 	using eigen = SEQ;
// };

template <typename...T>
struct TypeBox
{
	static constexpr auto size = sizeof...(T);
	using indices = n_seq_t<sizeof...(T)>;

	template <size_t N>
	using type = getType<N, T...>;
};

template <typename RET, bool REDUCE, typename ... Ops>
struct COMPUTATION : _COMPUTATION
{
	//14.07 this is not necessary any longer
	static_assert(((Ops::rep) | ... | REP::NONE) != REP::BOTH,
				  "Single computation cannot contain operators requiring different spaces (X and P).");


	using returnType = std::remove_const_t<std::remove_reference_t<RET>>;
	using types = std::remove_const_t<COMPUTATION>;

	static_assert(std::is_convertible_v<returnType, double> || sizeof(returnType) % sizeof(double) == 0, "The size of returnType is not a multiple of double.");

	static constexpr bool reduce = REDUCE;

	static constexpr size_t returnTypeSize = std::is_convertible_v<returnType, double> ? 1 : (sizeof(returnType) / sizeof(double));

	//Contains custom value in case of ALLOCS_EIGENSTATES base class
	using type_seq = n_seq_t<sizeof...(Ops)>;
	static constexpr size_t returnCount = (sizeof...(Ops));
	static constexpr size_t sizeInBuffer = returnCount * returnTypeSize;

	template <typename Op>
	static constexpr bool holdsOp = (std::is_same_v<Op, Ops> || ... || false);
	template <typename Op>
	static constexpr size_t bufferIndexOp = Index_v<Op, Ops...>*returnTypeSize;

	// using bufferOffsets = accumulate<mult_seq_t <returnTypeSize, n_seq_t<returnCount>>>;
	using bufferOffsets = mult_seq_t <returnTypeSize, n_seq_t<returnCount>>;
	static constexpr bool late = (Ops::late || ... || false);
	static constexpr REP rep = (Ops::rep | ... | REP::NONE);
	static constexpr inline std::string_view formatting = FMT_DOUBLE;
	static constexpr inline std::string_view format = repeat_v<returnCount, formatting>;


	// template <typename WHEN> static constexpr bool canRun = true;
	// template <MODE M, REP R, OPTIMS opts> inline static void prepare() {}
	// template <MODE M, REP R, OPTIMS opts> inline static void forerunner() {}
};

template <class Op> struct AVG : COMPUTATION <double, true, Op>
{
	static constexpr std::string_view name = "AVG_";
};

template <class ... Op> struct FLUX : COMPUTATION <double, true, Op...>
{
	static constexpr std::string_view name = "FLX_";
};

using ZOA_FLUX_3D = FLUX<N2S, N2D, N2T, S2D, S2T, S2CAP, D2CAP, T2CAP>;
using ZOA_FLUX_2D = FLUX<N2S, N2D, S2D, S2CAP, D2CAP>;


template <uind dir, class Op> struct DERIVATIVE : COMPUTATION <double, true, Op>
{
	static constexpr std::string_view name = "DER_";
};

template <typename ... Args> struct AUXILLARY_VALUES :COMPUTATION<double, false, Args...>
{
	static constexpr std::string_view name = "";
};


struct Sum
{
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = true;
};
struct Change
{
	static constexpr REP rep = REP::NONE;
	static constexpr bool late = true;
};


template <class ... Op>
struct SUM : COMPUTATION<std::common_type_t< typename Op::returnType...>, (Op::reduce || ...), Sum >
{

};

template <class Op> struct CHANGE : COMPUTATION<typename Op::returnType, false, Change>
{
};

template <typename ... Args> struct VALUE : COMPUTATION<double, false, Args...>
{
	static constexpr std::string_view name = "";
};

template <typename ... Args> struct OPERATION : COMPUTATION<double, false>
{
};

/*
template <typename...Ts>
struct Dumps
{
	static constexpr auto size = sizeof...(Ts);
	using type = Dumps;

	template <REP R, typename WHEN>
	static inline void run()
	{
		(runEach<Ts, R, WHEN>(), ...);
	}

	template <typename T, REP R, typename WHEN>
	static inline void runEach()
	{

	}
}; */


/*
template <size_t pos, class RetT = COMP, typename  typename COMP::returnT>
	inline void storeInBuffer(RetT val)
	{
		// logInfo("would store %s at %td", typeid(COMP).name(), pos);
		constexpr bool usingReduceBuffer = usesReduceBuffer<M, COMP>;
		// logInfo("About to stack... %td %d %g", pos, usingReduceBuffer, val);

		if constexpr (std::is_same_v<std::remove_const_t<RetT>, double> || std::is_convertible_v<RetT, double>)
		{
			// logInfo("here %d", xbufferCurrentLine);
			if constexpr (usingReduceBuffer) rbuffer[rbufferCurrentLine + pos] = val;
			else if (!MPI::pID)
			{
				// logInfo("Stacked at %td ", xbufferCurrentLine + pos);
				xbuffer[xbufferCurrentLine + pos] = val;
			}
			// logInfo("here2");
		}
		else if constexpr (std::is_same_v<std::remove_const_t<RetT>, cxd>)
		{
			if constexpr (usingReduceBuffer)
			{
				rbuffer[rbufferCurrentLine + pos] = val.real();
				rbuffer[rbufferCurrentLine + pos + 1] = val.imag();
			}
			else if (!MPI::pID)
			{
				xbuffer[xbufferCurrentLine + pos] = val.real();
				xbuffer[xbufferCurrentLine + pos + 1] = val.imag();
			}
		}
		else
		{
			constexpr size_t RetTsize = sizeof(RetT) / sizeof(double);
			for (int i = 0; i < RetTsize; i++)
			{
				if constexpr (usingReduceBuffer)rbuffer[rbufferCurrentLine + pos + i] = val[i];
				else if (!MPI::pID)  xbuffer[xbufferCurrentLine + pos + i] = val[i];
			}
		}
	} */
// template <typename...T>
// struct Run :TypeBox<T...>
// {

// 	using base = TypeBox<T...>;
// 	using indices = typename base::indices;
// 	// static constexpr indices {};
// 	static constexpr MODE modes[] = { T::mode... };
// 	// static constexpr auto runnableIndices = filter_seq<runnableQ, n_seq<sizeof...(T)>>{};

// 	// template <MODE M>
// 	// constexpr static bool canRunMode = ((size_t(M) & MODES) || MODES == 0);

// private:
// 	static inline void cancel(size_t I)
// 	{
// 		logWarning("Routine %td will not be executed, because mode %s is disabled via compilation option -MODES", I, modeName(modes[I]));
// 	}

// 	template <size_t...I>
// 	static inline void run(int argc, char* argv[], seq<I...>)
// 	{
// 		([&]
// 		 {
// 			//  if constexpr (canRunMode<T::mode>)
// 			//  {
// 			 T routine{ };
// 			 routine.template run<I>(argc, argv);
// 		//  }
// 		//  else cancel(I);
// 		 }(), ...);
// 	}
// public:
// 	static inline void run(int argc, char* argv[])
// 	{
// 		run(argc, argv, indices{});
// 	}
// };

// ForEach(comps, [&](auto index) {
			// 	using comp_type = tuple_element_t<index, decltype(comps)>;
			// 	constexpr auto comp = get<index>(comps);
			// 	if constexpr (0 < tuple_size_v < decltype(comp.types) >)
			// 	{
			// 		ForEach([&](auto comp.types,  subIndex) -> void
			// 				{

			// 					logQuick(comp.formatting(),
			// 							 (usesReduceBuffer<M, comp_type> ?
			// 							  rbuffer[rbufferLastLine + subIndex + getPositionInBuffer<RI, comp_type>()] :
			// 							  xbuffer[xbufferLastLine + subIndex + getPositionInBuffer<RI, comp_type>()]));

			// 				});
			// 	}});

			// constexpr auto aux_tup = filter_by_base<_AUXILLARY_VALUES>(getTasks<RI>());
			// ForEach(aux_tup, [&](auto index) -> void
			// 		{
			// 			constexpr auto aux = get<index>(aux_tup);
			// 			ForEach(aux.types, [&](auto subIndex) -> void
			// 					{
			// 						get<subIndex>(aux.types).log();
			// 					});
			// 		});

