#include "myutils.h"
#include "ppcaligner.h"
#include "abcxyz.h"

void PPCAligner::SetQuery(PDBChain &Q)
	{
	Q.CheckPPCMotifCoords();
	m_Q = &Q;
	}

void PPCAligner::SetRef(const PDBChain &R)
	{
	m_R = &R;
	asserta(m_R->m_MotifPosVec.size() == 3);
	}

double PPCAligner::GetMotifRMSD()
	{
	uint QueryPosA = m_Q->m_MotifPosVec[A];
	uint QueryPosB = m_Q->m_MotifPosVec[B];
	uint QueryPosC = m_Q->m_MotifPosVec[C];

	uint RefPosA = m_R->m_MotifPosVec[A];
	uint RefPosB = m_R->m_MotifPosVec[B];
	uint RefPosC = m_R->m_MotifPosVec[C];

	double RMSD2A = GetRMSD2Segment(QueryPosA, RefPosA, AL);
	double RMSD2B = GetRMSD2Segment(QueryPosB, RefPosB, BL);
	double RMSD2C = GetRMSD2Segment(QueryPosC, RefPosC, CL);
	double RMSD2 = RMSD2A + RMSD2B + RMSD2C;
	double Mean = (RMSD2A + RMSD2B + RMSD2C) / (AL + BL + CL);
	double RMSD = sqrt(Mean);
	return RMSD;
	}

double PPCAligner::GetRMSD2Segment(uint QPos, uint RPos, uint n)
	{
	double Sum = 0;
	for (uint i = 0; i < n; ++i)
		{
		m_Q->GetPt(QPos + i, m_QPt);
		m_R->GetPt(RPos + i, m_RPt);
		double d2 = GetDist2(m_QPt, m_RPt);
		Sum += d2;
		}
	return Sum;
	}

void PPCAligner::Align(TSHit &Hit) const
	{
	asserta(m_Q != 0 && m_R != 0);
	Hit.m_Query = m_Q;
	Hit.m_Ref = m_R;
	}
