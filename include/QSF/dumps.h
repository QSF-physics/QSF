
struct _DUMP {};

template <typename WHAT, typename WHEN>
struct DUMP :_DUMP
{
	using type = WHAT;
	WHEN when;
	DIMS dims;
	REP rep;
	double maxCoord;
	constexpr DUMP() {};
	constexpr DUMP(WHEN when, DIMS dims, REP rep, double maxCoord) :
		when(when), dims(dims), rep(rep), maxCoord(maxCoord) {}

	template <ind RI, REP R> inline void dump() const;
};

struct _PSI {};
struct _POT {};

template <REP R, DIMS D, size_t B = 1, size_t MC = 100000>
struct DUMP_PSI
{
	REP rep = R;
	DIMS dims = D;
	size_t bin = B;
	size_t maxCoord = MC;
	// template <ind RI, REP R> inline void dump() const;
};

template <REP R, DIMS D, size_t B = 1, size_t MC = 100000>
struct DUMP_POT
{
	REP rep = R;
	DIMS dims = D;
	size_t bin = B;
	size_t maxCoord = MC;
	// template <ind RI, REP R> inline void dump() const;
};
