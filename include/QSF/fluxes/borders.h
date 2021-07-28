#include "borders_g.h"
#include "borders_t.h"

/** STANDARD TWO AND THREE ELECTRON BORDERS
 * 			(Y direction)
 * 		0 		--->		  n
 *      L3	  L2 L1  R1 R2	  R3
 * 		|	  |  |	 |	|	  |
 * 				  A2
 * 0 	 _____________________  <-- p_CAP 		X
 *  	|		 |2S |		  |
 * |	|	2D	 |___|	1D	  | <-- p_NS		d
 * |	|________|	 |________| <-- p_SD		i
 * | A3 | 3S  |	   N	| 1S  |	A1 				r
 * V	|⎻⎻⎻⎻⎻⎻⎻⎻|   |⎻⎻⎻⎻⎻⎻⎻⎻| <-- p_SD
 * 		|	3D	 |⎻⎻⎻|	4D	  | <-- p_NS
 * 		|		 |4S |		  |
 * n 	 ⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻⎻  <-- p_CAP	n
 * 		^	  ^	 ^A4 ^	^     ^
 * 		|p_CAP|	 |	 |	|	  | p_CAP
 * 			  |	 |	 |	| p_NS
 * 			  |	 |	 | p_SD
 *            |  | p_SD
 * 			  | p_NS 					**/


// static ind p_NS = 10;
// static ind p_SD, p_DT, p_CAP = 10;
// p_NS: position of borders between N and S regions (standard: 12.5)
static ind p_NS = (ind)(12.5 / 0.5);
// p_SD: position of borders between S and D regions (standard: 7)
static ind p_SD = (ind)(7. / 0.5);
// p_DT: position of borders between D and T regions (standard: 5)
static ind p_DT = (ind)(5. / 0.5);
// p_CAP = (ind)(n2 - 1.0 * PROP::CAPlength / 0.5);// CAP_nodes;
static ind p_CAP = (ind)(20.0 / 0.5);

// Simplest flux border, takes size [au] 
template <int size>
struct BOX :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view size_str = num_to_string<size>;
	static constexpr std::string_view base_name = "BOX";
	static constexpr std::string_view name = join_v<base_name, num_to_string<size>>;

	static inline borVec<1> border(ind i) //1e
	{
		// double len = size / dx;
		return borVec<1>(PointB2(i, size));
	}

	static inline borVec<2> border(ind i, ind j)
	{
		// double size = size / dx;
		return borVec<2>(
			SegmentBi2j2(j, i, size, 0, size),
			SegmentBi2j2(i, j, size, 0, size)
			);
	}

	static inline borVec<3> border(ind i, ind j, ind k)
	{
		// double size = size / dx;
		return borVec<3>(
			PatchBi2j2k2(j, k, i, size, 0, size, 0, size),
			PatchBi2j2k2(k, i, j, size, 0, size, 0, size),
			PatchBi2j2k2(i, j, k, size, 0, size, 0, size));
	}
};

//1e, 2e and 3e borders
struct S2CAP :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "S2CAP";
	static inline borVec<2> border(ind i, ind j) //2e
	{
		return borVec<2>(SegmentBi2j2(j, i, p_CAP, 0, p_SD),
						 SegmentBi2j2(i, j, p_CAP, 0, p_SD));
	}
	static inline borVec<3> border(ind i, ind j, ind k) //3e
	{
		return borVec<3>(
			PatchBi2j2k2(j, k, i, p_CAP, 0, p_DT, 0, p_SD) + PatchBi2j2k2(k, j, i, p_CAP, 0, p_DT, 0, p_SD),
			PatchBi2j2k2(k, i, j, p_CAP, 0, p_DT, 0, p_SD) + PatchBi2j2k2(i, k, j, p_CAP, 0, p_DT, 0, p_SD),
			PatchBi2j2k2(i, j, k, p_CAP, 0, p_DT, 0, p_SD) + PatchBi2j2k2(j, i, k, p_CAP, 0, p_DT, 0, p_SD));
	}
};
//2e and 3e borders
struct D2CAP :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "D2CAP";
	static inline borVec<2> border(ind i, ind j)
	{
		return borVec<2>(SegmentBi2j2(j, i, p_CAP, p_SD, p_CAP),
						 SegmentBi2j2(i, j, p_CAP, p_SD, p_CAP));
	}
	static inline borVec<3> border(ind i, ind j, ind k)
	{
		return borVec<3>(
			PatchBi2j2k2(j, k, i, p_CAP, p_SD, p_CAP, 0, p_DT) + PatchBi2j2k2(k, j, i, p_CAP, p_SD, p_CAP, 0, p_DT),
			PatchBi2j2k2(k, i, j, p_CAP, p_SD, p_CAP, 0, p_DT) + PatchBi2j2k2(i, k, j, p_CAP, p_SD, p_CAP, 0, p_DT),
			PatchBi2j2k2(i, j, k, p_CAP, p_SD, p_CAP, 0, p_DT) + PatchBi2j2k2(j, i, k, p_CAP, p_SD, p_CAP, 0, p_DT));
	}
};
struct N2S :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "N2S";
	static inline borVec<2> border(ind i, ind j) //2e
	{
		return borVec<2>(SegmentBi2j2(j, i, p_NS, 0, p_SD),
						 SegmentBi2j2(i, j, p_NS, 0, p_SD));
	}
	static inline borVec<3> border(ind i, ind j, ind k) //3e
	{
		return borVec<3>(PatchBi2j2k2(j, k, i, p_NS, 0, p_SD, 0, p_DT - 1) + PatchBi2j2k2(k, j, i, p_NS, 0, p_SD, 0, p_DT - 1),
						 PatchBi2j2k2(k, i, j, p_NS, 0, p_SD, 0, p_DT - 1) + PatchBi2j2k2(i, k, j, p_NS, 0, p_SD, 0, p_DT - 1),
						 PatchBi2j2k2(i, j, k, p_NS, 0, p_SD, 0, p_DT - 1) + PatchBi2j2k2(j, i, k, p_NS, 0, p_SD, 0, p_DT - 1));
	}
};
struct N2D :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "N2D";
	static inline borVec<2> border(ind i, ind j) //2e
	{
		return borVec<2>(SegmentBi2j2(j, i, p_SD, p_SD, p_NS),
						 SegmentBi2j2(i, j, p_SD, p_SD, p_NS));
	}
	static inline borVec<3> border(ind i, ind j, ind k) //3e
	{
		return borVec<3>(PatchBi2j2k2(j, k, i, p_SD, p_SD, p_NS, 0, p_DT - 1) + PatchBi2j2k2(k, j, i, p_SD, p_SD, p_NS, 0, p_DT - 1),
						 PatchBi2j2k2(k, i, j, p_SD, p_SD, p_NS, 0, p_DT - 1) + PatchBi2j2k2(i, k, j, p_SD, p_SD, p_NS, 0, p_DT - 1),
						 PatchBi2j2k2(i, j, k, p_SD, p_SD, p_NS, 0, p_DT - 1) + PatchBi2j2k2(j, i, k, p_SD, p_SD, p_NS, 0, p_DT - 1));
	}
};

struct S2D :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "S2D";
	static inline borVec<2> border(ind i, ind j) //2e
	{
		return borVec<2>(SegmentBi2j2(j, i, p_SD, p_NS, p_CAP),
						 SegmentBi2j2(i, j, p_SD, p_NS, p_CAP));
	}
	static inline borVec<3> border(ind i, ind j, ind k) {
		return borVec<3>(PatchBi2j2k2(j, k, i, p_SD, p_NS, p_CAP, 0, p_DT) + PatchBi2j2k2(k, j, i, p_SD, p_NS, p_CAP, 0, p_DT),
						 PatchBi2j2k2(k, i, j, p_SD, p_NS, p_CAP, 0, p_DT) + PatchBi2j2k2(i, k, j, p_SD, p_NS, p_CAP, 0, p_DT),
						 PatchBi2j2k2(i, j, k, p_SD, p_NS, p_CAP, 0, p_DT) + PatchBi2j2k2(j, i, k, p_SD, p_NS, p_CAP, 0, p_DT));
	}
};

// 3e only
struct N2T :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "N2T";
	static inline borVec<3> border(ind i, ind j, ind k)
	{
		return borVec<3>(PatchBi2j2k2(j, k, i, p_DT, p_DT - 1, p_NS, p_DT - 1, p_SD) + PatchBi2j2k2(k, j, i, p_DT, p_DT - 1, p_NS, p_DT - 1, p_SD),
						 PatchBi2j2k2(k, i, j, p_DT, p_DT - 1, p_NS, p_DT - 1, p_SD) + PatchBi2j2k2(i, k, j, p_DT, p_DT - 1, p_NS, p_DT - 1, p_SD),
						 PatchBi2j2k2(i, j, k, p_DT, p_DT - 1, p_NS, p_DT - 1, p_SD) + PatchBi2j2k2(j, i, k, p_DT, p_DT - 1, p_NS, p_DT - 1, p_SD));
	}
};
struct S2T :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "S2T";
	static inline borVec<3> border(ind i, ind j, ind k)
	{
		return borVec<3>(PatchBi2j2k2(j, k, i, p_DT, p_NS, p_CAP, p_DT, p_SD) + PatchBi2j2k2(k, j, i, p_DT, p_NS, p_CAP, p_DT, p_SD),
						 PatchBi2j2k2(k, i, j, p_DT, p_NS, p_CAP, p_DT, p_SD) + PatchBi2j2k2(i, k, j, p_DT, p_NS, p_CAP, p_DT, p_SD),
						 PatchBi2j2k2(i, j, k, p_DT, p_NS, p_CAP, p_DT, p_SD) + PatchBi2j2k2(j, i, k, p_DT, p_NS, p_CAP, p_DT, p_SD));
	}
};
struct D2T :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "D2T";
	static inline borVec<3> border(ind i, ind j, ind k)
	{
		return borVec<3>(PatchBi2j2k2(j, k, i, p_DT, p_SD, p_CAP, p_SD, p_CAP),
						 PatchBi2j2k2(k, i, j, p_DT, p_SD, p_CAP, p_SD, p_CAP),
						 PatchBi2j2k2(i, j, k, p_DT, p_SD, p_CAP, p_SD, p_CAP));
	}
};
struct T2CAP :_Operator
{
	static constexpr REP rep = REP::X;
	static constexpr bool late = false;
	static constexpr std::string_view name = "T2CAP";
	static inline borVec<3> border(ind i, ind j, ind k)
	{
		return borVec<3>(
			PatchBi2j2k2(j, k, i, p_CAP, p_DT, p_CAP, p_DT, p_CAP),
			PatchBi2j2k2(k, i, j, p_CAP, p_DT, p_CAP, p_DT, p_CAP),
			PatchBi2j2k2(i, j, k, p_CAP, p_DT, p_CAP, p_DT, p_CAP));
	}
};
