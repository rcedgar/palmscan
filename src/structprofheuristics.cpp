#include "myutils.h"
#include "structprof.h"
#include "abcxyz.h"

/***
D1 = max dist to G in 50aa window after motif C
  (start of pre-D helix).
D2 = min dist to G in 50aa window after D1.
D3 = min dist to A in +/- 10aa window around D2.
***/
uint StructProf::FindMofifD_Hueuristics() const
	{
	const uint L = m_Chain->GetSeqLength();

	uint Pos_aD = m_Chain->GetMotifPos(A) + 3;
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	uint Pos_cD = m_Chain->GetMotifPos(C) + 3;

	double Dist1 = DBL_MAX;
	uint Lo1 = Pos_cD;
	uint Hi1 = Lo1 + 50;
	if (Hi1 >= L)
		Hi1 = L - 1;
	uint D1 = SearchDist(Pos_bG, Lo1, Hi1, true, 5.0, Dist1);

	double Dist2 = DBL_MAX;
	uint Lo2 = D1;
	uint Hi2 = D1 + 50;
	if (Hi2 >= L)
		Hi2 = L - 1;
	uint D2 = SearchDist(Pos_bG, Lo2, Hi2, false, 5.0, Dist2, true);
	if (D2 < 3)
		return UINT_MAX;

	uint CN = GetCavityNumber(D2);
	if (CN > 6)
		return UINT_MAX;

	return D2;
	}

uint StructProf::FindMofifE_Hueuristics(uint Pos_MotifD) const
	{
	if (Pos_MotifD == UINT_MAX)
		return UINT_MAX;
	const uint L = m_Chain->GetSeqLength();
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	uint Lo1 = Pos_MotifD + 8;
	uint Hi1 = Lo1 + 25;
	if (Hi1 >= L)
		Hi1 = L - 1;
	double Dist1 = DBL_MAX;
	uint Pos1 = SearchDist(Pos_bG, Lo1, Hi1, true, 5.0, Dist1);
	if (Pos1 >= Hi1 - 5)
		return UINT_MAX;
	double DistE;
	uint PosE = SearchDist(Pos_bG, Pos1+1, Hi1, false, 5.0, DistE);
	if (PosE < 3)
		return UINT_MAX;
	return PosE - 3;
	}

// PosF2 + 2 = 'R'
uint StructProf::FindMofifF2_Hueuristics(uint Pos_MotifA) const
	{
	if (Pos_MotifA == UINT_MAX)
		return UINT_MAX;
	if (Pos_MotifA < 25)
		return UINT_MAX;
	uint Pos_bG = m_Chain->GetMotifPos(B) + 1;
	double Dist = DBL_MAX;
	uint Lo1 = (Pos_MotifA < 100 ? 0 : Pos_MotifA - 100);
	uint Hi1 = Pos_MotifA - 10;
	uint Pos1 = SearchDist(Pos_bG, Lo1, Hi1, false, 99.0, Dist);

	int iPos;
	const char *Seq = m_Chain->m_Seq.c_str();

#define CheckOffset(Offset)	\
	{ \
	iPos = int(Pos1) + Offset; \
	if (iPos >= 0 && Seq[iPos + 2] == 'R') \
		return uint(iPos); \
	}

	CheckOffset(-6)
	CheckOffset(-5)
	CheckOffset(-4)

#undef CheckOffset

	if (Pos1 >= 5)
		return Pos1 - 5;

	return UINT_MAX;
	}
