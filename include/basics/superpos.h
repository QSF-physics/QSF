
template <size_t... I> struct superpos : seq<I...>
{
	double weight[sizeof...(I)];
	static constexpr auto pos = n_seq<sizeof...(I)>;
	template<typename ... Args>
	constexpr superpos(Args... args) :seq<I...>(), weight{ args... }{}
};