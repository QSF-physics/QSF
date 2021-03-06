#pragma once

/* ::::::::::::::::::: Base Implementation ::::::::::::::::::: */
template <REP first, REP... rest>
struct _SplitBase
{
	// double mults[sizeof...(rest) + 1];
	static constexpr REP firstREP = first;
	//Those below need to be converted back to REP
	using restREP = seq<size_t(rest)...>;
	using fullREP = seq<size_t(first), size_t(rest)...>;
	using indices = up_to_t<sizeof...(rest)>;
	static constexpr size_t restCount = sizeof...(rest);
	static constexpr size_t baseSplitsNumber = 1 + restCount;
	static constexpr bool canJoin = (first == (rest, ...));
};

template <REP first, REP middle, REP last>
struct Split3Base :_SplitBase<first, middle, last>
{
	static constexpr double mults[3] = { 0.5,1.0,0.5 };
	static constexpr std::string_view name = "BASE";
	using base = _SplitBase<first, middle, last>;
	using base::restCount;
	using base::baseSplitsNumber;
	using base::canJoin;
};

/* ::::::::::::::::::: Chain Implementation ::::::::::::::::::: */
template <typename Base, typename ... Coeffs>
struct _SplitChain
{
	static constexpr size_t N = sizeof...(Coeffs);
	static constexpr std::string_view times = "x";
	static constexpr std::string_view pl = "_";
	static constexpr std::string_view N_str = num_to_string<N>;
	static constexpr std::string_view name = join_v<pl, Base::name, times, N_str>;


	using reps = std::conditional_t < Base::canJoin, concat_seq < seq<size_t(Base::firstREP)>, repeat_seq<N, typename Base::restREP>>, repeat_seq<N, typename Base::fullREP>>;
	using splits = n_seq_t<reps::size>;

	static constexpr auto pack() noexcept
	{
		// double arr[reps::size];
		std::array<double, reps::size> arr{};
		for (uind i = 0; i < reps::size; i++) arr[i] = 0;
		int i = 0;
		([&i, &arr]
		 {
			 for (uind j = 0; j < Base::baseSplitsNumber; j++)
				 arr[i + j] += Base::mults[j] * Coeffs::value;
			 i += (Base::baseSplitsNumber - (Base::canJoin ? 1 : 0));
		 }(), ...);
		return arr;

	}
	static constexpr auto mults = pack();
	template <uind N>
	static constexpr REP rep = REP(nth_seq_elem<N>(reps{}));
};

template <typename Base, size_t Order>
struct MultiProductSplit
{
	static_assert(Order > 0, "Order of MultiProductSplit method must be 1>=0.");
	static constexpr REP firstREP = Base::firstREP;
	static constexpr std::string_view name = "MP";
	using ChainExpander = from_to_t<1, Order>;

	template <size_t N, size_t Ns>
	struct SplitCoeff
	{
		static constexpr double value = 1.0 / N;
	};

	//this is not necessairly equal to the number of real splits
	template <size_t N>
	using SplitExpander = from_to_t<1, N>;

	template <size_t I, class = SplitExpander<I>, class = ChainExpander >
	struct Chain;
	template <size_t I, size_t ... Js, size_t ... Is>
	struct Chain < I, seq<Js...>, seq<Is...>> : _SplitChain<Base, SplitCoeff<I, Js>...>
	{
		static constexpr double value = (...* (Is != I ? double(I * I) / (double(I * I) - double(Is * Is)) : 1.0));
	};
	// static constexpr size_t ChainCount = Order;
	// template <typename N>
	// using chain = Chain<N>;
};