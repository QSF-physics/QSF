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

template <typename RET, bool red, typename ... Ops>
struct COMPUTATION : _COMPUTATION
{
	//14.07 this is not necessary any longer
	static_assert(((Ops::rep) | ... | REP::NONE) != REP::BOTH,
				  "Single computation cannot contain operators requiring different spaces (X and P).");


	using returnType = std::remove_const_t<std::remove_reference_t<RET>>;
	using types = std::remove_const_t<COMPUTATION>;

	static_assert(std::is_convertible_v<returnType, double> || sizeof(returnType) % sizeof(double) == 0, "The size of returnType is not a multiple of double. Hint: Try casting your type to ");

	static constexpr bool reduce = red;

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

//Not stored in buffers, but shown during logging. Useful for ETAOperator
struct _AUXILLARY_VALUES {};
template <typename ... Args> struct AUXILLARY_VALUES :_AUXILLARY_VALUES
{
	std::tuple<Args...> types;
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


template <class ... Op> struct SUM :COMPUTATION<std::common_type_t< typename Op::returnType...>, (Op::reduce || ...), Sum >
{
	// std::tuple<Op...> types;
};

template <class Op> struct CHANGE : COMPUTATION<typename Op::returnType, Op::reduce, Change>
{
	// std::tuple<Op> types;
};

template <typename ... Args> struct VALUE : COMPUTATION<double, false, Args...>
{
	static constexpr std::string_view name = "";
};

template <typename ... Args> struct OPERATION : COMPUTATION<double, false>
{
};

struct BufferedOutputsBase
{
	int comp_interval;			// How often time dependent info will be writen
	int log_interval;			// How often logs will be written to stdout 
	FILE* file_dat = nullptr; 	// for writing the results of time evolution
};

template <bool binary, typename...Ts>
struct BufferedOutputs : BufferedOutputsBase, TypeBox<Ts...>
{
	using TypeBox<Ts...>::size;

	template <typename T>
	static constexpr bool usesReduceBuffer = T::reduce;

	template < bool usingReduceBuffer>
	static constexpr size_t sizeInBuffer = (((usesReduceBuffer< Ts> == usingReduceBuffer) ? Ts::sizeInBuffer : 0) + ... + 0);

	template < typename T>
	static constexpr size_t offset()
	{
		bool found = false;
		size_t index = 0;
		([&] { if (!found) {
			if constexpr (std::is_same_v<Ts, T>) found = true;
			else if (usesReduceBuffer<Ts> == usesReduceBuffer<T>)
				index += Ts::sizeInBuffer;
		}}(), ...);
		return index;
	}

	template < typename T>
	using pos = offset_seq_t<offset<T>(), typename T::bufferOffsets>;


	// using all_pos = concat_all_seq<pos<Ts>...>;

	uind bufferHeight;

	double* rbuffer;
	uind rbufferSize;
	uind rbufferCurrentLine;
	uind rbufferLastLine;

	double* xbuffer;
	uind xbufferSize;
	uind xbufferCurrentLine;
	uind xbufferLastLine;

	template < typename T>
	double* record()
	{
		if constexpr (usesReduceBuffer<T>)
			return rbuffer + rbufferCurrentLine;
		else return xbuffer + xbufferCurrentLine;
	}
	template <MODE M>
	void init(uind PASS, std::string_view name)
	{
		logInfo("comp_interval %d", comp_interval);
		logInfo("log_interval %d", log_interval);

		// (printf("size in buffers: %d %td %s", Ts::sizeInBuffer, TypeBox<Ts...>::size, typeid(Ts).name()), ...);
		bufferHeight = log_interval / comp_interval;
		rbufferSize = sizeInBuffer<true> *bufferHeight;
		xbufferSize = sizeInBuffer<false> *bufferHeight;

		rbufferLastLine = (bufferHeight - 1) * sizeInBuffer<true>;
		xbufferLastLine = (bufferHeight - 1) * sizeInBuffer<false>;
		rbufferCurrentLine = rbufferLastLine;
		xbufferCurrentLine = xbufferLastLine;
		logBUFFER("REDUCABLE BUFFER SIZE: %tdx%d, NORMAL BUFFER SIZE: %tdx%d", sizeInBuffer<true>, bufferHeight, sizeInBuffer<false>, bufferHeight);
		if (rbufferSize > 0) rbuffer = new double[rbufferSize];
		if (!MPI::pID)
		{
			if (xbufferSize > 0) xbuffer = new double[xbufferSize];
		}

		if (comp_interval > 0)
			file_dat = openOut<IO_ATTR::WRITE >(name, PASS, binary);
	}

	BufferedOutputs(Section& settings)
	{
		inipp::get_value(settings, "comp_interval", comp_interval);
		inipp::get_value(settings, "log_interval", log_interval);
	}

	BufferedOutputs(BufferedOutputsBase bob) :BufferedOutputsBase(bob)
	{}

	~BufferedOutputs()
	{
		if (rbufferSize > 0) { logSETUP("Destroying rbuffer"); delete[] rbuffer; }
		if (!MPI::pID) if (xbufferSize > 0) { logSETUP("Destroying xbuffer"); delete[] xbuffer; }
		closeFile(file_dat);
	}

	template < class COMP>
	auto getLastValue()
	{
		// return 11.0;
		return usesReduceBuffer<COMP>
			? rbuffer[rbufferCurrentLine + offset<COMP>()]
			: xbuffer[xbufferCurrentLine + offset<COMP>()];
	}


	template < class COMP, class...Op>
	void store(Op...val)
	{
		using returnType = typename COMP::returnType;
		constexpr bool usingReduceBuffer = usesReduceBuffer<COMP>;
		//If the computation *WOULD* use reduce buffer in different mode reduce it immediataly


		size_t pos = offset<COMP>();

		// logInfo("About to stack %s with val %g at pos %td using %s", typeid(COMP).name(), val..., pos, usingReduceBuffer ? "rbuff" : "xbuff");

		([&] {

			if constexpr (std::is_same_v<returnType, double> || std::is_convertible_v<returnType, double>)
			{
				// logInfo("vals %g %g", val...);
				if constexpr (usingReduceBuffer) rbuffer[rbufferCurrentLine + pos] = val;
				else if (!MPI::pID)
				{
					// logInfo("Stacked at %td ", xbufferCurrentLine + pos);
					xbuffer[xbufferCurrentLine + pos] = val;
				}
				// logInfo("here2");
			}
			else if constexpr (std::is_same_v<returnType, cxd>)
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
				constexpr size_t RetTsize = sizeof(returnType) / sizeof(double);
				for (int i = 0; i < RetTsize; i++)
				{
					if constexpr (usingReduceBuffer)rbuffer[rbufferCurrentLine + pos + i] = val[i];
					else if (!MPI::pID)  xbuffer[xbufferCurrentLine + pos + i] = val[i];
				}
			}
			pos += COMP::returnTypeSize;
		 }(), ...);
	}
	// template <typename WHEN, REP R, class PROP>
	// inline void compute(const PROP& propagator)
	// {
	// 	((Ts::template canRun<R, WHEN> ?
	// 	  computeEach<PROP, WHEN, R, Ts, typename Ts::returnT>
	// 	  (propagator, Ts{}, (typename Ts::types) {}, pos<Ts>{}) : ditch()), ...);
	// }




	// template <class PROP, typename WHEN, REP R, typename T, typename retT, typename... COMP, size_t...Is>
	// inline void computeEach(const PROP& propagator, T&& comp, COMPUTATION<retT, Op...>&&, seq<	// Is...>&&)
	// {
	   // 	// T::template forerunner<R, opt>();
	   // 	// logInfo("compute:");
	   // 	// Timings::measure::start(comp.name);
	   // 	(storeInBuffer < Is, usesReduceBuffer<T>, retT>(
	   // 		propagator->template calc<R, Op>()), ...);
	   // 		// comp.template calc<R, opt, Op>()), ...);
	   // 		// Timings::measure::stop(comp.name);

	   // 	// runEach<R, opt, >(T{});
	   // }

	static inline void setupComputations()
	{
		// if (optims & PRECOMP_COORDS) precomputeCoords();

		// ForEach(comps, [&](auto index)
		// 		{
		// 			// constexpr auto comps = getComputations<RI, _COMPUTATION>();
		// 			constexpr auto x = get<index>(comps);

		// 			if constexpr (x.template canRun<getPropagator<RI>.startsIn, EARLY>())
		// 				x.template prepare<mode, REP::X, optims>();

		// 			if constexpr (x.template canRun<(REP::BOTH ^ getPropagator<RI>.startsIn), LATE>())
		// 				x.template prepare<mode, REP::P, optims>();

		// 		});
	}


	template <typename ...T>
	void logQuick(const std::string_view format, T... data)
	{
		LOG_INLINE(format.data(), data...);
	}


	void writeCaptions()
	{

		// auto names = logNames<RI>();
		// logInfo("%s", names.data());
		// constexpr auto comps = getComputations<RI>();
		// if constexpr (tuple_size_v < decltype(comps) > > 0)
		// {
			// LOG_INLINE_START(__LOG_NC);

			// ForEach(comps, [](auto index)
		// 			{
		// 				constexpr auto comps = getComputations<RI>();
		// 				// using type = tuple_element_t<index, decltype(comps)>;
		// 				constexpr auto item = get<index>(comps);
		// 				// using item_type = decltype(item);
		// 				if constexpr (0 < tuple_size_v < decltype(item.types) >)
		// 				{
		// 					ForEach(item.types, [&](auto subIndex) -> void
		// 							{
		// 								constexpr auto val = get<subIndex>(item.types);
		// 								if constexpr (is_base_of_v<_Operator, decltype(val)>)
		// 									logQuick("|%13s ", join_v<item.name, val.name>);
		// 								if constexpr (!is_base_of_v<_Operator, decltype(val)>)
		// 									LOG_INLINE("|%11s%2zu ", item.name.data(), val());
		// 									//TODO: should be get<subIndex>(item.types)()

		// 							});
		// 				}
		// 			});

		// 	constexpr auto aux_tup = filter_by_base<_AUXILLARY_VALUES>(getTasks<RI>());
		// 	ForEach(aux_tup, [&](auto index) -> void
		// 			{
		// 				constexpr auto aux = get<index>(aux_tup);
		// 				ForEach(aux.types, [&](auto subIndex) -> void
		// 						{
		// 							LOG_INLINE("|%13s ", get<subIndex>(aux.types).name.data());
		// 						});
		// 			});
		// 	LOG_INLINE_END();
	// } // }
	}



	template <typename T, size_t ... I>
	inline void log(seq <I...>)
	{
		// logInfo("logging %g from %d",((usesReduceBuffer<T> ? rbuffer + rbufferLastLine : xbuffer + xbufferLastLine) + I))
		LOG_INLINE(T::format.data(), *((usesReduceBuffer<T> ? (rbuffer + rbufferLastLine) : xbuffer + xbufferLastLine) + I)...);
	}
	// template <typename T>
	// inline void log(seq<I...>)

	// template <typename T>
	// inline void log()
	// {
	// 	// logInfo(format.data(), data...);
	// }



	template <typename XBUFF, typename RBUFF>
	void writeDataBinaryHeader(int columns, XBUFF xBuff, RBUFF rBuff)
	{
		// if (file_dat != nullptr)
		// {
		// 	int dim = DIM;
		// 	fwrite(&dim, sizeof(int), 1, file_dat);
		// 	fwrite(&n, sizeof(int), 1, file_dat);
		// 	// int mode = M;
		// 	// fwrite(&mode, sizeof(int), 1, file_dat);
		// 	fwrite(&columns, sizeof(int), 1, file_dat);

		// 	ForEach(xBuff, [&](auto sub)
		// 			{
		// 				dispatchComp(file_dat, get<sub>(xBuff));
		// 			});
		// 	ForEach(rBuff, [&](auto sub)
		// 			{
		// 				dispatchComp(file_dat, get<sub>(rBuff));
		// 			});
		// }
	}

	template <MODE M, WHEN when>
	inline void logAll()
	{
		if constexpr (when == WHEN::AT_END)
		{
			rbufferLastLine = rbufferCurrentLine;
			xbufferLastLine = xbufferCurrentLine;
		}
		// constexpr auto comps = getComputations<RI>();
		if constexpr (bool(size))
		{
			LOG_INLINE_START(__LOG_NC);
			(log<Ts>(pos< Ts>{}), ...);

			LOG_INLINE_END();
		}

		if (!MPI::pID) //WRITE DAT FILE
		{
			int i;
			int end;
			if constexpr (binary && when == WHEN::AT_START)
			{
				// constexpr auto x_comp = getComputations<RI, _XBUFFER>();
				// constexpr auto r_comp = getComputations<RI, _RBUFFER>();
				// writeDataBinaryHeader(sizeInBuffer<false> +sizeInBuffer<true>, x_comp, r_comp);
			}
			if constexpr (when == WHEN::AT_START)
			{
				i = bufferHeight - 1;
				end = bufferHeight;
			}
			else if constexpr (when == WHEN::AT_END)
			{
				i = 0;
				end = (xbufferCurrentLine / sizeInBuffer<false>) + 1;
			}
			else
			{
				i = 0;
				end = bufferHeight;
			}
			if constexpr (binary)
			{
				for (; i < end; i++)
				{
					fwrite(xbuffer + sizeInBuffer<false> *i, sizeof(double), sizeInBuffer<false>, file_dat);
					fwrite(rbuffer + sizeInBuffer<true> *i, sizeof(double), sizeInBuffer<true>, file_dat);
				}
			}
			else for (; i < end; i++)
			{
				for (int j = 0; j < sizeInBuffer<false>; j++)
					fprintf(file_dat, "%15.5g ", xbuffer[i * sizeInBuffer<false> +j]);
				for (int j = 0; j < sizeInBuffer<true>; j++)
					fprintf(file_dat, "%15.5g ", rbuffer[i * sizeInBuffer<true> +j]);
				fprintf(file_dat, FMT_END);
			}
			fflush(file_dat);
		}
	}

	inline void reduce()
	{
		if (rbufferSize > 0)
		{
			if (!MPI::pID) MPI_Reduce(MPI_IN_PLACE, rbuffer, rbufferSize,
									  MPI_DOUBLE, MPI_SUM, 0, MPI::rComm);
			else MPI_Reduce(rbuffer, rbuffer, rbufferSize,
							MPI_DOUBLE, MPI_SUM, 0, MPI::rComm);
		}
	}

	inline void reduceLine()
	{
		if (rbufferSize > 0)
		{
			if (!MPI::pID) MPI_Reduce(MPI_IN_PLACE, rbuffer + rbufferCurrentLine, sizeInBuffer<true>, MPI_DOUBLE, MPI_SUM, 0, MPI::rComm);
			else MPI_Reduce(rbuffer + rbufferCurrentLine, rbuffer + rbufferCurrentLine, sizeInBuffer<true>, MPI_DOUBLE, MPI_SUM, 0, MPI::rComm);
		}
	}

	template <MODE M, WHEN when>
	inline void logOrPass(ind step)
	{

		if ((log_interval > 0 && (step % log_interval == 0)) || when == WHEN::AT_END)
		{
			if (M == MODE::IM) reduceLine();
			if (M == MODE::RE) reduce();

			if (!MPI::pID) logAll<M, when>();
			rbufferCurrentLine = 0;
			xbufferCurrentLine = 0;
		}
		else if (comp_interval > 0 && step % comp_interval == 0)
		{
			if (M == MODE::IM) reduceLine();
			rbufferCurrentLine += sizeInBuffer<true>;
			xbufferCurrentLine += sizeInBuffer<false>;
		}
	}

	void test()
	{
		logTestFatal(log_interval >= comp_interval, "log_interval (%d) should be greater than comp_interval (%d)", log_interval, comp_interval);
		logTestFatal(log_interval % comp_interval == 0, "log_interval (%d) should be divisable by comp_interval (%d)", log_interval, comp_interval);
	}
};

template <typename ... Ts>
using BufferedBinaryOutputs = BufferedOutputs<true, Ts...>;
template <typename ... Ts>
using BufferedTextOutputs = BufferedOutputs<false, Ts...>;


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

