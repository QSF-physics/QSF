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
	using base::canJoin;
	using base::baseSplitsNumber;
	using base::restCount;
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
		for (int i = 0; i < reps::size; i++) arr[i] = 0;
		int i = 0;
		([&i, &arr]
		 {
			 for (int j = 0; j < Base::baseSplitsNumber; j++)
				 arr[i + j] += Base::mults[j] * Coeffs::value;
			 i += (Base::baseSplitsNumber - (Base::canJoin ? 1 : 0));
		 }(), ...);
		return arr;

	}
	static constexpr auto mults = pack();
};

template <typename Base, size_t Order>
struct MultiProductSplit
{
	static constexpr std::string_view name = "MP";
	using ChainExpander = from_to_t<1, Order>;

	template <size_t N>
	struct SplitCoeff
	{
		static constexpr double value = 1.0 / N;
	};

	//this is not necessairly equal to the number of real splits
	template <size_t N>
	using InternalSplitExp = ChainExpander;

	template <size_t I, class = InternalSplitExp<I> >
	struct Chain;
	template <size_t I, size_t ... Is>
	struct Chain < I, seq<Is...>> : _SplitChain<Base, SplitCoeff<Is>...>
	{
		static constexpr double value = (...* (Is != I ? double(I * I) / (double(I * I) - double(Is * Is)) : 1.0));
	};
	// static constexpr size_t ChainCount = Order;
	// template <typename N>
	// using chain = Chain<N>;
};



// template <typename Base, size_t Order, template <class, size_t> class Type>
// struct OperatorSplit
// {
// 	using SplitBase = Base;
// 	using SplitType = Type<Base, Order>;

// 	static constexpr size_t numberOfChains = SplitType::numberOfChains;

// 	static constexpr size_t numberOfChains = SplitType::numberOfChains;
// 	template <size_t N>
// 	using chain = typename SplitType::template chain<N>;

// 	template <size N>
// 	static constexpr chainSize = chain<N>::size;

// 	static constexpr double mults[chain_type::size];
// 	using SumExpander = n_seq_t<count>;
// 	static constexpr REP startsIn = Base::firstREP;
// 	static constexpr std::string_view name = join_v<Chains::name...>;
// 	size_t sizes[sizeof...(Chains)];
// };

// template <size_t ... Is>
// constexpr double Vandermondian(size_t i, seq<Is...>)
// {
// 	return (...* (Is != i ? double(i * i) / (double(i * i) - double(Is * Is)) : 1.0));
// }

// _SplitSum<Base, ChainCoeff<Base, from_to_t<1, Is>>...>
// template <typename Base, typename expand>
// struct MultiProductChain {};

// template <typename Base, auto ... Is>
// struct MultiProductChain < Base, seq<Is...>> : _SplitChain<sizeof...(Is), Base>
// {
// 	using Super = _SplitChain<sizeof...(Is), Base>;
// 	using expander = typename Super::expander;
// 	static constexpr double baseMult = 1.0 / Super::bases;
// 	double coeff;
// 	constexpr MultiProductChain(double coeff) :coeff(coeff),
// 		Super(baseMult, 10.0, expander{}) {}
// 		// Super::splits{ typename Base::first_split_t{baseMult,5.0}, tuple_element_t<Is,typename Base::splits_t>{baseMult, 5.0}... } {}
// };

// template <typename Base, size_t Order, typename expand = from_to_t<1, Order>> struct _MultiproductSum {};

// template <typename Base, size_t Order, size_t...Is>
// struct _MultiproductSum <Base, Order, seq<Is...>> : _SplitSum<Base, MultiProductChainCoeff<Base, from_to_t<1, Is>>...>
// {
// 	using Super = _SplitSum<Base, MultiProductChain<Base, from_to_t<1, Is>>...>;
// 	static constexpr size_t Order = sizeof...(Is);
// 	constexpr _MultiproductSum() :
// 		Super{ MultiProductChain<Base,  from_to_t<1, Is>>{ Vandermondian(Is, from_to<1,Order>) }... } {}
// };
// namespace detail
// {
// 	template <typename Base, size_t Order, typename expand = from_to_t<1, Order>>
// 	using MultiproductSum = _MultiproductSum<Base, expand>;

// }

// template <typename Base, size_t Order = 1>
// using MultiproductSum = detail::MultiproductSum<Base, 1>;