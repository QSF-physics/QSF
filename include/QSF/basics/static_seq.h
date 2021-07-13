
//index in sequence parameter pack
template <size_t T, size_t H, size_t... I>
constexpr size_t index_in_seq()
{
	if constexpr (bool(sizeof ...(I))) return T == H ? 0 : 1 + index_in_seq<T, I...>();
	else return 0;
}

//basic non-negative integer (index) seq type
template<size_t... Ns> struct seq
{
	using type = seq;
	static constexpr size_t size = sizeof...(Ns);
};

//concat two seq types
template <typename SeqA, typename SeqB> struct concat;
template <size_t... N1, size_t... N2> struct concat<seq<N1...>, seq<N2...>>
{
	using type = seq<N1..., N2...>;
};
template <typename SeqA, typename SeqB>
using concat_seq = typename concat<SeqA, SeqB>::type;
template <typename H, typename... T>
struct concat_all : concat_seq<H, typename concat_all<T...>::type> {};

//concat any number of seq types
template<auto... Ms> struct concat_all<seq<Ms...>> { using type = seq<Ms...>; };
template <typename ...Seqs> using concat_all_seq = typename concat_all<Seqs...>::type;

template <size_t ... Is>
constexpr bool is_in(size_t i, seq<Is...>)
{
	return ((i == Is) || ...);
}
template <typename SEQ, size_t ... I>
constexpr bool contains_seq(SEQ sq, seq<I...>)
{
	return (is_in(I, sq) && ...);
}

//UNIQ SEQUENCE
namespace uniq_impl
{
	template<class, auto... Ns> struct uniq;

	template<auto... Ms, auto N, auto... Ns> struct uniq<seq<Ms...>, N, Ns...> :
		std::conditional_t<(... || (N == Ms)),
		uniq<seq<Ms...>, Ns...>,
		uniq<seq<Ms..., N>, Ns...>> {};

	template<auto... Ms> struct uniq<seq<Ms...>>
	{
		using type = seq<Ms...>;
	};

	template <class> struct uniq_seq;
	template <size_t... I> struct uniq_seq<seq<I...>> : uniq<seq<>, I...> {};
};
template<class SEQ> using uniq_seq = typename uniq_impl::uniq_seq<SEQ>::type;

template<size_t... Ns> using uniq = typename uniq_impl::uniq<seq<>, Ns...>::type;

namespace seq_gen_impl
{
	template <class, class> struct merge;

	template <size_t... I1, size_t... I2>
	struct merge<seq<I1...>, seq<I2...>> : seq<I1..., (sizeof...(I1) + I2)...> {};

	template <size_t N> struct seq_gen : merge<typename seq_gen<N / 2>::type, typename seq_gen<N - N / 2>::type> {};

	template<> struct seq_gen<0> : seq<> { };
	template<> struct seq_gen<1> : seq<0> { };

	template <size_t O, class> struct seq_offset;
	template <size_t O, size_t ... I> struct seq_offset<O, seq<I...>> : seq<(O + I)...> {};

	template <size_t M, class> struct seq_mult;
	template <size_t M, size_t ... I> struct seq_mult<M, seq<I...>> : seq<(M* I)...> {};

	template <size_t S, size_t E, typename = std::enable_if_t<(S <= E)>> struct seq_gen_from : seq_offset<S, typename seq_gen<E - S>::type> { static_assert(S < E, "S<E"); };
};

//FILTER SEQUENCE
namespace filter_impl
{
	// template <template <typename> typename Pred, auto... Ns>
	// struct filter
	// {
	// 	using type
	// }
	// template<class, auto... Ns> struct filter;
	template <template <auto> typename Pred, class> struct filter;

	template<template <auto> typename Pred, auto... Ns>
	struct filter <Pred, seq<Ns...>> :
		concat_all_seq < std::conditional_t<Pred<Ns>::value, seq<Ns>, seq<>>...> {};

	// template <template <auto> typename Pred, size_t... I>
	// struct filter_seq<seq<I...>> : filter<Pred, I...> {};
	// template<auto... Ms> struct filter<seq<Ms...>>
	// {
	// 	using type = seq<Ms...>;
	// };

	// template <class> struct filter_seq;
	// template <size_t... I> struct filter_seq<seq<I...>> : filter<seq<>, I...> {};
	// template <auto val, class Pred>
	// constexpr auto filter_single(seq<val>, Pred pred)
	// {
		// if constexpr (pred(val)) return seq<val>;
		// else return seq<>
	// }
	// template<class Pred, size_t...Ns>
	// constexpr auto filter(seq<Ns...>, Pred pred)
	// {
	// 	if constexpr (sizeof...(Ns)>0)
	// 		return concat_all_seq<
	// }
};
template<template <auto> typename Pred, class SEQ> using filter_seq = typename filter_impl::filter<Pred, SEQ>::type;
// template<size_t... Ns> using filter = typename uniq_impl::filter<seq<>, Ns...>::type;


//TYPES
template <size_t N> using n_seq_t = typename seq_gen_impl::seq_gen<N>::type;
template <size_t N> using up_to_t = typename seq_gen_impl::seq_gen<N + 1>::type;
template <size_t S, size_t E> using from_to_t = typename seq_gen_impl::seq_gen_from<S, E + 1>::type;

template <size_t M, typename SEQ> using mult_seq_t = typename seq_gen_impl::seq_mult<M, SEQ>::type;
template <size_t O, typename SEQ> using offset_seq_t = typename seq_gen_impl::seq_offset<O, SEQ>::type;
//OBJS
template <size_t N> constexpr static auto n_seq = n_seq_t<N>{};
template <size_t N> constexpr static auto up_to = up_to_t<N>{};
template <size_t S, size_t E> constexpr static auto from_to = from_to_t<S, E>{};


namespace impl
{
	template <class, size_t, class>
	struct accumulate;

	template <size_t...Is, size_t S, size_t H, size_t...Js>
	struct accumulate <seq<Is...>, S, seq<H, Js...>> : accumulate<seq <Is..., S + H>, S + H, seq<Js...>> {};

	template <size_t ...Is, size_t S>
	struct accumulate < seq<Is...>, S, seq<>>
	{
		using type = seq<Is...>;
	};
}
// template<class>
// struct accumulate;

template<typename SEQ>
using accumulate = typename impl::accumulate<seq<>, 0, SEQ>::type;

namespace detail
{
	template <size_t N, typename SEQ>
	struct repeat_seq
	{
		template <typename = n_seq_t<N>> struct impl;

		template <size_t... Is> struct impl<seq<Is...>>
		{
			template <size_t>
			using wrap = SEQ;

			using type = concat_all_seq <wrap<Is>... >;
		};
		using type = typename impl<>::type;
	};
}
template <size_t N, typename SEQ>
using repeat_seq = typename detail::repeat_seq<N, SEQ>::type;

/* For types */
template <typename T, typename... Ts>
struct Index;

template <typename T>
struct Index<T> : std::integral_constant<size_t, 0> {};

template <typename T, typename... Ts>
struct Index<T, T, Ts...> : std::integral_constant<size_t, 0> {};

template <typename T, typename U, typename... Ts>
struct Index<T, U, Ts...> : std::integral_constant<size_t, 1 + Index<T, Ts...>::value> {};

template <typename T, typename... Ts>
constexpr size_t Index_v = Index<T, Ts...>::value;