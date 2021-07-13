
// template <size_t... I> struct superpos : seq<I...>
// {
// 	double weight[sizeof...(I)];
// 	static constexpr auto pos = up_to<sizeof...(I) - 1>;
// 	template<typename ... Args>
// 	constexpr superpos(Args... args) :seq<I...>(), weight{ args... }{};
// };

template <bool customSource, typename SP = superpos<0>> struct source
{
	static constexpr bool custom = customSource;
	size_t index = 0;
	SP superpositions = superpos<0>{ 1.0 };
	constexpr source() = default;
	// constexpr source(size_t index): index(index){};
	constexpr source(size_t index, SP superpositions = superpos<0>{ 1.0 }) :superpositions(superpositions), index(index) {}
	constexpr source(SP superpositions, size_t index = 0) : superpositions(superpositions), index(index) {}
};
// template <typename SP> source(SP) -> source<false, SP>;
template <typename SP, class = typename std::enable_if<!std::is_integral<SP>::value>::type> source(SP)->source<false, SP>;
template <typename SP> source(size_t, SP)->source<true, SP>;
template <typename SP> source(SP, size_t)->source<true, SP >;
source(size_t)->source<true>;
source()->source<false>;

