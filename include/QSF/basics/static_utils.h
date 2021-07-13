#pragma once
#include <tuple>

auto is_even = [](int _in) {return _in % 2 == 0;};

template <typename...>
inline constexpr auto is_unique = std::true_type{};

template <typename T, typename... Rest>
inline constexpr auto is_unique<T, Rest...> = std::bool_constant<
	(!std::is_same_v<T, Rest> && ...) && is_unique<Rest...>
>{};


// template <typename ...Args>
// using areCopyConstructible = typename std::conjunction<std::is_copy_constructible<Args>...>::type;

namespace detail
{

	template<size_t, class = void, class...>
	struct getType;

	template<size_t N, class T, class... Rest>
	struct getType
	{
		static_assert(!is_same_v<T, void>, "Can't get n-th element from empty container");
		using type = typename getType<N - 1, Rest...>::type;
	};

	template<class T, class... Rest>
	struct getType<0, T, Rest...>
	{
		using type = T;
	};
}

template<size_t N, typename... Args>
using getType = typename detail::getType<N, Args...>::type;

template <template <class> class, template <class...> class, class...>
struct filter;

template <template <class> class Pred, template <class...> class Variadic>
struct filter<Pred, Variadic>
{
	using type = Variadic<>;
};

template <template <class> class Pred,
	template <class...> class Variadic,
	class T, class... Ts>
	struct filter<Pred, Variadic, T, Ts...>
{
	template <class, class> struct Cons;
	template <class Head, class... Tail>
	struct Cons<Head, Variadic<Tail...> >
	{
		using type = Variadic<Head, Tail...>;
	};

	using type = std::conditional_t<
		Pred<T>::value,
		typename Cons<T, typename filter<Pred, Variadic, Ts...>::type>::type,
		typename filter<Pred, Variadic, Ts...>::type >;
};


template <class T, template <class...> class Template>
struct pred : std::false_type {};
template <template <class...> class Template, class... Args>
struct pred<Template<Args...>, Template> : std::true_type {};

// template <template <typename> typename Pred, typename... Ts> struct mfilter_all<Pred, tuple<Ts...>>
template <template <typename> typename Pred, typename... Ts> struct mfilter_all
{
	using type = decltype(tuple_cat(
		std::conditional_t<Pred<Ts>::value, std::tuple<Ts>, std::tuple<>>{}...
	));
};



// int main() {
// 	static_assert(std::is_same<
// 				  filter<std::is_integral, std::tuple, int, float, long>::type,
// 				  std::tuple<int, long> >::value, "");
// }

			// template <template <typename...> typename Template, typename Tup>
			// constexpr auto filter_tuple(Tup&& tup)
			// {
			// 	return std::apply([](auto...ts) {
			// 		return std::tuple_cat(std::conditional_t<(is_spec<decltype(ts), Template>::value),
			// 							  std::tuple<decltype(ts)>,
			// 							  std::tuple<>>{}...);
			// 					  }, tup);
			// }


template <class Tuple, class F, std::size_t... I>
constexpr F for_each_impl(Tuple&& t, F&& f, std::index_sequence<I...>)
{
	return (void)std::initializer_list<int>{(std::forward<F>(f)(std::get<I>(std::forward<Tuple>(t))), 0)...}, f;
}

template <class Tuple, class F>
constexpr F for_each(Tuple&& t, F&& f)
{
	return for_each_impl(std::forward<Tuple>(t), std::forward<F>(f),
						 std::make_index_sequence<std::tuple_size_v<std::remove_reference_t<Tuple>>>{});
}

template<typename Tuple, typename Predicate>
constexpr size_t find_if(Tuple&& tuple, Predicate pred)
{
	size_t index = std::tuple_size<std::remove_reference_t<Tuple>>::value;
	size_t currentIndex = 0;
	bool found = false;
	for_each(tuple, [&](auto&& value)
			 {
				 if (!found && pred(value))
				 {
					 index = currentIndex;
					 found = true;
				 }
				 ++currentIndex;
			 });
	return index;
}



///EFFICIENT FOR EACH https://stackoverflow.com/questions/56937863/constexpr-lambda-argument
// Unluckily fails with gcc-8.3

template <class Func, std::size_t... index>
static inline constexpr void ForEach(Func&& f, std::index_sequence<index...>)
{
	(f(std::integral_constant<std::size_t, index>()), ...);
}


template <class Func, typename ... Tp>
static inline constexpr void ForEach(std::tuple<Tp...>, Func&& f)
{
	if constexpr ((sizeof ... (Tp)) > 0)
		ForEach(f, std::make_index_sequence<sizeof ... (Tp)>());
}



// template<size_t I = 0, typename F, typename... Tp>
// constexpr void it_comp_names(F f, tuple<Tp...> t)
// {
// 	for (size_t j = 0; j < get<I>(t).returnCount;j++) f(get<I>(t).name());
// 	if constexpr (I + 1 != sizeof...(Tp)) it_comp_names<I + 1>(f, t);
// }
// template<size_t I = 0, typename F, typename... Tp>
// constexpr void it_comp_formats(F f, tuple<Tp...> t)
// {
// 	for (size_t j = 0; j < get<I>(t).returnCount;j++) f(get<I>(t).formatting());
// 	if constexpr (I + 1 != sizeof...(Tp)) it_comp_formats<I + 1>(f, t);
// }



//general type concat
template <typename T, typename ...>
struct cat
{
	using type = T;
};

template <typename ... Ts1, typename ... Ts2, typename ... Ts3>
struct cat<std::tuple<Ts1...>, std::tuple<Ts2...>, Ts3...>
	: public cat<std::tuple<Ts1..., Ts2...>, Ts3...>
{ };

template <typename ...Ts >
using cat_t = typename cat <Ts...>::type;

// template <size_t N, typename Tup>
// class NTimesTupleTypes
// {
// 	template <typename = n_seq_t<N>> struct impl;

// 	template <size_t... Is> struct impl<seq<Is...>>
// 	{
// 		template <size_t>
// 		using wrap = Tup;

// 		using type = cat_t <wrap<Is>... >;
// 	};
// public:
// 	using type = typename impl<>::type;
// 	static constexpr auto size = std::tuple_size_v<type>;
// };

// template<size_t I = 0, typename... Tp>
// void print2(std::tuple<Tp...> t) {
// 	cout << get<I>(t) << " ";
// 	// do things
// 	if constexpr (I + 1 != sizeof...(Tp))
// 		print2<I + 1>(t);
// }


// template <class T>
// struct usesReduceBuffer : std::false_type {};

	// template <class T>
	// struct rBufferT : std::conjunction<is_base_of_v<_MPI_REDUCE, T>, !is_base_of_v<_MPI_IMMEDIATE, T>> {};
	// template <class T>
	// struct xBufferT : std::disjunction<!is_base_of_v<_MPI_REDUCE, T>, is_base_of_v<_MPI_IMMEDIATE, T>> {};



		// template<typename T, typename C, std::size_t I>
		// struct tuple_index_r;

		// template<typename H, typename ...R, typename C, std::size_t I>
		// struct tuple_index_r<std::tuple<H, R...>, C, I>
		// 	: public std::conditional<std::is_same<C, H>::value,
		// 	std::integral_constant<std::size_t, I>,
		// 	tuple_index_r<std::tuple<R...>, C, I + 1>>::type
		// {};

		// template<typename C, std::size_t I>
		// struct tuple_index_r<std::tuple<>, C, I>
		// {};

		// template<typename T, typename C>
		// struct tuple_index
		// 	: public std::integral_constant<std::size_t, tuple_index_r<T, C, 0>::value> {};

		//SECOND TRY
		// template <class T, class Tuple>
		// struct tuple_index
		// {
		// 	static_assert(!is_same_v<Tuple, tuple<>>, "Could not find `T` in given `Tuple`");
		// };

		// template <class T, class... Ts>
		// struct tuple_index<tuple<T, Ts...>, T> {
		// 	static const size_t value = 0;
		// };

		// template <class T, class U, class... Ts>
		// struct tuple_index<tuple<U, Ts...>, T> {
		// 	static const size_t value = 1 + tuple_index<tuple<Ts...>, T>::value;
		// };

		// template <class... Args>
		// struct type_list
		// {
		//    template <std::size_t N>
		//    using type = typename std::tuple_element<N, std::tuple<Args...>>::type;
		// };

		// int main()
		// {
		//    std::cout << std::boolalpha;
		//    type_list<int, char, bool>::type<2> x = true;
		//    std::cout << x << '\n';
		// }





// // template <IND RI, typename TU, typename T>
// // constexpr auto getLastRoutine(const TU& ROUTS, const T& value)
// // {
// // 	for (size_t i = RI; i > 0; i--)
// // 		if (_andQ<RI>(ROUTS, value))
// // 			return get<i>(ROUTS);
// // 	return get<RI>(ROUTS);
// // }

// template<size_t I = 0, typename FuncT, typename... Tp>
// inline typename enable_if<I == sizeof...(Tp), void>::type
// for_each2(tuple<Tp...>, FuncT) // Unused arguments are given no names.
// { }

// template<size_t I = 0, typename FuncT, typename... Tp>
// inline typename enable_if < I < sizeof...(Tp), void>::type
// 	for_each2(tuple<Tp...> t, FuncT f)
// {
// 	f(get<I>(t));
// 	for_each2<I + 1, FuncT, Tp...>(t, f);
// }



// template<size_t I, typename T, typename Pred, std::enable_if_t<I == 0, int> = 0>
// constexpr T get_element_helper(T const* arr, Pred&& pred) {
// 	return pred(arr[0]) ? arr[0] : get_element_helper<0>(arr + 1, pred);
// }
// template<size_t I, typename T, typename Pred, std::enable_if_t<I != 0, int> = 0>
// constexpr T get_element_helper(T const* arr, Pred&& pred) {
// 	return pred(arr[0]) ? get_element_helper<I - 1>(arr + 1, pred) : get_element_helper<I>(arr + 1, pred);
// }

// template<typename T, size_t N, typename Pred, size_t ...Is>
// constexpr std::array<T, sizeof...(Is)> filter_array_helper(std::array<T, N> const& arr, Pred&& pred, std::index_sequence<Is...>*) {
// 	return { get_element_helper<Is>(arr.data(), pred)... };
// }
// template<size_t M, typename T, size_t N, typename Pred>
// constexpr std::array<T, M> filter_array(std::array<T, N> const& arr, Pred&& pred) {
// 	return filter_array_helper(arr, std::forward<Pred>(pred), (std::make_index_sequence<M>*)nullptr);
// }

// template<size_t I = 0, typename F, typename... Tp>
// constexpr auto _get_from_tuple(tuple<Tp...> t, F pred) {
// 	if constexpr (I + 1 != sizeof...(Tp) && !F(get<I>(t)))
// 		return _get_from_tuple<I + 1>(t, pred);
// 	else return get<I>(t); //error
// }

// template<size_t I = 0, typename T, typename... Tp>
// constexpr auto get_from_tuple(tuple<Tp...> t) {
// 	if constexpr (I + 1 != sizeof...(Tp) && !is_base_of_v<T, decltype(get<I>(t))>)
// 		return get_from_tuple<I + 1>(t);
// 	else return get<I>(t); //error
// }

// template<size_t I = 0, typename... Tp, typename T>
// void tuple_elems_based_on(tuple<Tp...> t) 
// {
// 	if constexpr (I + 1 != sizeof...(Tp)) 
// 	{
// 		return tuple_cat((get<I>(t))...);

// 	}
// 	else return tuple<>;
// }
// std::conditional_t<(is_same_v<T, decltype(ts)>),
// 							  std::tuple<decltype(ts)>,
// 							  std::tuple<>>{}...

// template<size_t I = 0, typename F, typename... Tp>
// constexpr size_t getIndexInTuple(tuple<Tp...> t, F pred) {
// 	if constexpr (I + 1 != sizeof...(Tp) && !F(get<I>(t)))
// 		return getIndexInTuple<I + 1>(t, pred);
// 	else return I;
// }
// template<size_t I = 0, typename... Tp>
// constexpr int get_size(tuple<Tp...> t)
// {
// 	if constexpr (I + 1 != sizeof...(Tp)) return 1 + get_size<I + 1>(t);
// }



// // func is enabled if all Ts... have the same type as T
// template<typename T, typename... Ts>
// std::enable_if_t<std::conjunction_v<std::is_same<T, Ts>...>>
// func(T, Ts...) {
//     std::cout << "all types in pack are T\n";
// }

// // otherwise
// template<typename T, typename... Ts>
// std::enable_if_t<!std::conjunction_v<std::is_same<T, Ts>...>>
// func(T, Ts...) {
//     std::cout << "not all types in pack are T\n";
// }

// int main() {
//     func(1, 2, 3);
//     func(1, 2, "hello!");
// }

// template<typename T, typename... Ts>
// constexpr bool contains()
// { return std::disjunction_v<std::is_same<T, Ts>...>; }

// static_assert(    contains<int,      bool, char, int, long>());
// static_assert(    contains<bool,     bool, char, int, long>());
// static_assert(    contains<long,     bool, char, int, long>());
// static_assert(not contains<unsigned, bool, char, int, long>());


// //C++17 does not have so we implement one easily
// template<typename T, size_t N, typename Pred>
// constexpr size_t count_if(array<T, N> const& arr, Pred&& pred)
// {
// 	auto siz = 0;
// 	for (const auto& v : arr) siz += (bool)pred(v);
// 	return siz;
// }


// template <class T, template <class...> class Template>
// struct is_spec : std::false_type {};
// template <template <class...> class Template, class... Args>
// struct is_spec<Template<Args...>, Template> : std::true_type {};


// template<typename, typename>
// struct tuple_holds;
// template<typename ...Ts, typename T>
// struct tuple_holds<std::tuple<Ts...>, T>
// 	: std::disjunction<std::is_same<Ts, T>...> {};

// // template<typename, typename> struct tuple_holds_child;
// // template<typename ...Ts, typename T> struct tuple_holds_child<std::tuple<Ts...>, T>
// 	// : std::disjunction<std::is_base_of<T, Ts>...> {};

// template <typename T, template <typename...> typename Template>
// struct tuple_holds_spec : std::false_type {};
// template<template <typename...> typename Template, typename ...Ts>
// struct tuple_holds_spec<std::tuple<Ts...>, Template>
// 	: std::disjunction<is_spec<remove_const_t<Ts>, Template>...> {};

// template< size_t I, typename Tuple_t, typename T>
// constexpr size_t tuple_index() {
// 	if constexpr (I < tuple_size<Tuple_t>::value)
// 	{
// 		using type = typename tuple_element<I, Tuple_t>::type;
// 		if constexpr (is_same_v<T, type>) return I;
// 		else return tuple_index<I + 1, Tuple_t, T >();
// 	}
// 	else return tuple_size<Tuple_t>::value;
// }

// template <auto Start, auto End, auto Inc, class F>
// constexpr void cfor(F&& f)
// {
// 	if constexpr (Start < End)
// 	{
// 		f(std::integral_constant<decltype(Start), Start>());
// 		cfor<Start + Inc, End, Inc>(f);
// 	}
// }
// 		// template <typename T, typename Tup>
// 		// constexpr auto filter_tuple(Tup&& tup)
// 		// {
// 		// 	return std::apply([](auto...ts) {
// 		// 		return std::tuple_cat(std::conditional_t<(is_same_v<T, decltype(ts)>),
// 		// 							  std::tuple<decltype(ts)>,
// 		// 							  std::tuple<>>{}...);
// 		// 					  }, tup);
// 		// }

// template <typename ... BASES, typename T>
// constexpr auto hasBase(const T t)
// {
// 	constexpr bool bo = (true && ... && is_base_of_v<BASES, T>);
// 	if constexpr (bo)  return tuple(t);
// 	else return tuple{};
// }

// template <typename ... WHEN, typename T>
// constexpr auto whenMatches(const T t)
// {
// 	constexpr bool bo = (false || ... || is_same_v<WHEN, decltype(t.when)>);
// 	if constexpr (bo)  return tuple(t);
// 	else return tuple{};
// }
// template <typename ... BASES, typename Tup>
// constexpr auto filter_by_base(const Tup tup)
// {
// 	return apply([](auto...ts) constexpr { return tuple_cat((hasBase<BASES...>(ts))...); }, tup);
// }
// template <typename... WHEN, typename Tup>
// constexpr auto filter_by_when(const Tup tup)
// {
// 	if constexpr (tuple_size_v<Tup> == 0)
// 		return tuple{};
// 	else
// 		return apply([](auto...ts) constexpr { return tuple_cat((whenMatches<WHEN...>(ts))...); }, tup);
// }
