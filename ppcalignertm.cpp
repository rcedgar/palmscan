#include "myutils.h"
#include "ppcaligner.h"
#include "abcxyz.h"

double AlignTmPpc(XDPMem &Mem, float **DPScoreMx,
  const PDBChain &Query, const PDBChain &Ref, string &Path);

//double PPCAligner::GetTMScore(const TSHit &Hit, string &Path) const
//	{
//	const PDBChain &Q = *Hit.m_Query;
//	const PDBChain &R = *Hit.m_Ref;
//	double TM = AlignTmPpc(Q, R, Path);
//	return TM;
//	}

double PPCAligner::GetTMScore()
	{
	asserta(m_Q != 0);
	asserta(m_R != 0);
	double TM = AlignTmPpc(m_Mem, m_DPScoreMx, *m_Q, *m_R, m_Path);
	return TM;
	}
