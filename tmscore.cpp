#include "myutils.h"
#include "trisearcher.h"
#include "abcxyz.h"

/***
	d0 = 1.24 cuberoot(L - 15) - 1.8  (eq.5)

	T_i = 1/{1 + (d_i/d0)^2}

	TM = 1/L Sum(i=1..L) T_i

When h = 0.75, eq. (4) can be well approximated by a 
simpler formula (eq. 5) which drops, for example,
from 6.4 to 2.3 Angstrom when LN changes from 300
to 50 residues.

tmscore_d0.py 300 = 6.36
tmscore_d0.py 50 = 2.26

Motif lengths 12+14+8 = 34

tmscore_d0.py 34 = 1.51
***/

static const double d0 = 1.51;
static const double d02 = d0*d0;

double TriSearcher::GetMotifsTM() const
	{
	uint QPosA, QPosB, QPosC;
	bool Found = GetTopHit(QPosA, QPosB, QPosC);
	if (!Found)
		return DBL_MAX;

	asserta(m_Ref != 0);
	asserta(m_Ref->m_MotifPosVec.size() == 3);

	uint RPosA = m_Ref->m_MotifPosVec[A];
	uint RPosB = m_Ref->m_MotifPosVec[B];
	uint RPosC = m_Ref->m_MotifPosVec[C];

	double Sum = 0;
	Sum += GetTMSum(QPosA, RPosA, AL);
	Sum += GetTMSum(QPosB, RPosB, BL);
	Sum += GetTMSum(QPosC, RPosC, CL);

	double TM = Sum/(AL + BL + CL);
	return TM;
	}

double TriSearcher::GetTMSum(uint QPos, uint RPos, uint n) const
	{
	vector<double> RefPt(3);
	vector<double> QPt(3);
	vector<double> QPtX(3);

	double Sum = 0;
	for (uint i = 0; i < n; ++i)
		{
		m_Query->GetPt(QPos + i, QPt);
		m_Ref->GetPt(RPos + i, RefPt);

		XFormPt(QPt, m_TriForm_t, m_TriForm_R, QPtX);

		double d2 = GetDist2(RefPt, QPtX);
		double T = 1.0/(1.0 + d2/d02);
		Sum += T;
		}
	return Sum;
	}
