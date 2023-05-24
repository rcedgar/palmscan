#include "myutils.h"
#include "trisearcher.h"
#include "abcxyz.h"

static const uint MAX_DP = 256;

#define TRACE	1

float ViterbiFastMem(XDPMem &Mem, float **ScoreMx,
  uint LA, uint LB, string &Path);

void TriSearcher::AllocDPScoreMx()
	{
	if (m_DPScoreMx != 0)
		return;

	m_DPScoreMx = myalloc(float *, MAX_DP);
	for (uint i = 0; i < MAX_DP; ++i)
		m_DPScoreMx[i] = myalloc(float, MAX_DP);
	}

void TriSearcher::AlignSeg(
  uint QPos, uint QLen, uint RPos, uint RLen,
  const vector<double> &t, const vector<vector<double> > &R,
  string &Path)
	{
	asserta(QLen > 0 && QLen < MAX_DP);
	asserta(RLen > 0 && RLen < MAX_DP);

	AllocDPScoreMx();

	vector<double> QPt(3);
	vector<double> QPtX(3);
	vector<double> RPt(3);
	for (uint i = 0; i < QLen; ++i)
		{
		m_Query->GetPt(QPos + i, QPt);
		XFormPt(QPt, t, R, QPtX);
		for (uint j = 0; j < RLen; ++j)
			{
			m_Ref->GetPt(RPos + j, RPt);

			float d = (float) GetDist(RPt, QPtX);
			m_DPScoreMx[i][j] = -d;
			}
		}
#if TRACE
	{
	Log("\n");
	Log("ScoreMx\n");
	Log("      ");
	for (uint j = 0; j < RLen; ++j)
		Log("   %8u", j);
	Log("\n");

	for (uint i = 0; i < QLen; ++i)
		{
		Log("[%4u]", i);
		for (uint j = 0; j < RLen; ++j)
			{
			float d = m_DPScoreMx[i][j];
			Log("  %9.3g", d);
			}
		Log("\n");
		}
	Log("\n");
	}
#endif
	ViterbiFastMem(m_XDPMem, m_DPScoreMx, QLen, RLen, Path);
#if TRACE
	Log("Path=%s\n", Path.c_str());
#endif
	const uint ColCount = SIZE(Path);
	uint i = QPos;
	uint j = RPos;
	const string &QSeq = m_Query->m_Seq;
	const string &RSeq = m_Ref->m_Seq;
	string QRow;
	string RRow;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c = Path[Col];
		if (c == 'M')
			{
			QRow += QSeq[i++];
			RRow += RSeq[j++];
			}
		else if (c == 'D')
			{
			QRow += QSeq[i++];
			RRow += '-';
			}
		else if (c == 'I')
			{
			QRow += '-';
			RRow += RSeq[j++];
			}
		}
#if TRACE
	{
	Log("QRow %s\n", QRow.c_str());
	Log("TRow %s\n", RRow.c_str());
	asserta(i == QPos + QLen);
	asserta(j == RPos + RLen);
	}
#endif
	}

bool TriSearcher::AlignPalm(uint QueryPosA, uint QueryPosB,
  uint QueryPosC, string &Path)
	{
	Path.clear();
	uint QPosA, QPosB, QPosC;
	bool Found = GetTopHit(QPosA, QPosB, QPosC);
	if (!Found)
		return false;

	asserta(m_Ref != 0);
	asserta(m_Ref->m_MotifPosVec.size() == 3);

	uint RPosA = m_Ref->m_MotifPosVec[A];
	uint RPosB = m_Ref->m_MotifPosVec[B];
	uint RPosC = m_Ref->m_MotifPosVec[C];

	uint V1StartQ = QPosA + AL;
	uint V1LenQ = QPosB - (QPosA + AL);

	uint V1StartR = RPosA + AL;
	uint V1LenR = RPosB - (RPosA + AL);

	string V1Path;
	AlignSeg(V1StartQ, V1LenQ, V1StartR, V1LenR,
	  m_TriForm_t, m_TriForm_R, V1Path);

	uint V2StartQ = QPosB + BL;
	uint V2LenQ = QPosC - (QPosB + BL);

	uint V2StartR = RPosB + BL;
	uint V2LenR = RPosC - (RPosB + BL);

	string V2Path;
	AlignSeg(V2StartQ, V2LenQ, V2StartR, V2LenR,
	  m_TriForm_t, m_TriForm_R, V2Path);

	int QL = m_Query->GetSeqLength();
	uint PalmLenQ = PDBChain::GetPalmPrintLength(QPosA, QPosC, QL);
	uint RL = SIZE(m_Ref->m_Seq);
	AlignSeg(QPosA, PalmLenQ, 0, RL,
	  m_TriForm_t, m_TriForm_R, Path);

	return true;
	}
