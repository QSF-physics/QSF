/*outwards oriented point boundary at iat*/
inline bool Rn(ind i, ind il, ind ir) { return (i >= il && i < ir); }

inline bool Rn2(ind i, ind il, ind ir) { return Rn(i, il, ir) || Rn(i, -ir, -il); }

/* iat - 1 because zero is not at n2*dx but at (n2-0.5)*dx */
inline double PointB(ind i, ind iat) { return i == iat - 1 ? Sign(i) : 0; }
inline double PointB2(ind i, ind iat) { return i == iat - 1 || i == -iat ? Sign(i) : 0; }
/*outwards oriented line boundary at jat from il to ir*/

inline double SegmentB(ind i, ind j, ind jat, ind il, ind ir)
{
	return Rn(i, il, ir) ? PointB(j, jat) : 0;
}
inline double SegmentBi2(ind i, ind j, ind jat, ind il, ind ir)
{
	return Rn2(i, il, ir) ? PointB(j, jat) : 0;
}
inline double SegmentBj2(ind i, ind j, ind jat, ind il, ind ir)
{
	return Rn(i, il, ir) ? PointB2(j, jat) : 0;
}
inline double SegmentBi2j2(ind i, ind j, ind jat, ind il, ind ir)
{
	return Rn2(i, il, ir) ? PointB2(j, jat) : 0;
}
/*outwards oriented plane boundary at kat from il to ir and jl to jr*/

inline double PatchB(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn(i, il, ir) ? SegmentB(j, k, kat, jl, jr) : 0;
}
inline double PatchBi2(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn2(i, il, ir) ? SegmentB(j, k, kat, jl, jr) : 0;
}
inline double PatchBj2(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn(i, il, ir) ? SegmentBi2(j, k, kat, jl, jr) : 0;
}
inline double PatchBk2(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn(i, il, ir) ? SegmentBj2(j, k, kat, jl, jr) : 0;
}
inline double PatchBi2j2(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn2(i, il, ir) ? SegmentBi2(j, k, kat, jl, jr) : 0;
}
inline double PatchBi2k2(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn2(i, il, ir) ? SegmentBj2(j, k, kat, jl, jr) : 0;
}
inline double PatchBj2k2(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn(i, il, ir) ? SegmentBi2j2(j, k, kat, jl, jr) : 0;
}
inline double PatchBi2j2k2(ind i, ind j, ind k, ind kat, ind il, ind ir, ind jl, ind jr)
{
	return Rn2(i, il, ir) ? SegmentBi2j2(j, k, kat, jl, jr) : 0;
}