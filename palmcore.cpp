#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

void PDBChain::GetPC(PDBChain &PC) const
	{
	PC.Clear();
	asserta(SIZE(m_MotifPosVec) == 3);

	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	uint L = GetSeqLength();
	asserta(PosB >= PosA + AL);
	asserta(PosC >= PosB + BL);
	asserta(PosC + CL <= L);

	uint Lo, Hi;
	GetPalmCoreCoords(PosA, PosC, L, Lo, Hi);

	PDBChain Trunc;
	TruncateChain(Lo, Hi, Trunc);
//	Trunc.GetTriFormChain_MotifCoords(PC);
	Trunc.GetTriFormChain_DGD(PC);

	uint PC_PosA = PosA - Lo;
	uint PC_PosB = PosB - Lo;
	uint PC_PosC = PosC - Lo;
	uint PCL = Hi - Lo + 1;

	PC.m_MotifPosVec.clear();
	PC.m_MotifPosVec.push_back(PC_PosA);
	PC.m_MotifPosVec.push_back(PC_PosB);
	PC.m_MotifPosVec.push_back(PC_PosC);

	string MotifA, MotifB, MotifC;
	GetSubSeq(PosA, AL, MotifA);
	GetSubSeq(PosB, BL, MotifB);
	GetSubSeq(PosC, CL, MotifC);

	PC.m_Label = m_Label;
	}
