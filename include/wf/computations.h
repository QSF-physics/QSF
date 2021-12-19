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
using ZOA_FLUX_2D_SIMPLE = FLUX<N2S, N2D, S2D, S2CAP, D2CAP>;
using ZOA_FLUX_2D = FLUX<N2S, N2D_SYM, N2D_ASYM, S2D_SYM, S2D_ASYM>;


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