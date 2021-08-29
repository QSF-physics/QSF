#include <array>
#include <string_view>
// // STRINGVIEW CONCAT
// // https ://stackoverflow.com/questions/38955940/how-to-concatenate-static-strings-at-compile-time
namespace detail
{
	using namespace std;

	template <string_view const&... Strs>
	struct join
	{
		// Join all strings into a single array of chars
		static constexpr auto impl() noexcept
		{
			constexpr size_t len = (Strs.size() + ... + 0);
			array<char, len + 1> arr{};
			auto append = [i = 0, &arr](auto const& s) mutable {
				for (auto c : s) arr[i++] = c;
			};
			(append(Strs), ...);
			arr[len] = 0;
			return arr;
		}
		// Give the joined string static storage
		static constexpr auto arr = impl();
		// View as a string_view
		static constexpr string_view value{ arr.data(), arr.size() - 1 };
	};

	template < size_t N, string_view const& Str>
	struct repeat
	{
		// Join all strings into a single array of chars
		static constexpr auto impl() noexcept
		{
			constexpr size_t len = Str.size() * N;
			array<char, len + 1> arr{};
			auto append = [i = 0, &arr](auto const& s) mutable {
				for (auto c : s) arr[i++] = c;
			};
			for (size_t j = 0; j < N; j++) append(Str);
			arr[len] = 0;
			return arr;
		}
		// Give the joined string static storage
		static constexpr auto arr = impl();
		// View as a string_view
		static constexpr string_view value{ arr.data(), arr.size() - 1 };
	};

	/* Nice simple int to string */
	template<size_t... digits>
	struct to_chars {
		static constexpr auto size = sizeof...(digits);
		static constexpr char value[] = { ('0' + digits)..., 0 };
	};

	template<size_t rem, size_t... digits>
	struct explode : explode<rem / 10, rem % 10, digits...> {};

	template<size_t... digits>
	struct explode<0, digits...> : to_chars<digits...> {};

	template <>
	struct explode<0> : to_chars<0> {};
}

// Helper to get the value out
template <std::string_view const&... Strs>
static constexpr auto join_v = detail::join<Strs...>::value;

// Helper to get the value out
template <size_t N, std::string_view const& Strs>
static constexpr auto repeat_v = detail::repeat<N, Strs>::value;

// constexpr auto num_to_string = detail::explode<num>::value;
template<size_t num>
static constexpr std::string_view num_to_string{ detail::explode<num>::value ,detail::explode<num>::size };