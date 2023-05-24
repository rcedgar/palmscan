#include "myutils.h"
#include "gspaligner.h"
#include "outputfiles.h"
#include "viterbi.h"

double GSPAligner::GetScore_Vecs(const vector<double> &QVec,
  const vector<double> &TVec) const
	{
	double Sum = 0;
	asserta(SIZE(QVec) == m_FeatureCount);
	asserta(SIZE(TVec) == m_FeatureCount);
	for (uint i = 0; i < m_FeatureCount; ++i)
		{
		double d = fabs(QVec[i] - TVec[i]);
		Sum += d;
		}
	double Score = 1 - Sum - m_Bias;
	return Score;
	}

double GSPAligner::GetScore_Pos(uint QPos, uint TPos) const
	{
	asserta(QPos < SIZE(m_Q->m_FeatureVec));
	asserta(TPos < SIZE(m_T->m_FeatureVec));
	const vector<double> &QVec = m_Q->m_FeatureVec[QPos];
	const vector<double> &TVec = m_T->m_FeatureVec[TPos];
	double Score = GetScore_Vecs(QVec, TVec);
	return Score;
	}

void GSPAligner::Alloc(uint QL, uint TL)
	{
	if (QL < m_MxSize && TL < m_MxSize)
		return;

	if (m_Mx != 0)
		{
		for (uint i = 0; i < m_MxSize; ++i)
			myfree(m_Mx[i]);
		myfree(m_Mx);
		}

	uint NewSize = max(QL, TL) + 100;
	m_Mx = myalloc(float *, NewSize);
	for (uint i = 0; i < NewSize; ++i)
		m_Mx[i] = myalloc(float, NewSize);
	}

void GSPAligner::MxToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	const uint QL = m_Q->GetSeqLength();
	const uint TL = m_T->GetSeqLength();
	for (uint QPos = 0; QPos < QL; ++QPos)
		{
		for (uint TPos = 0; TPos < TL; ++TPos)
			{
			double Score = m_Mx[QPos][TPos];
			fprintf(f, "%u\t%u\t%.4f\n",
			  QPos, TPos, Score);
			}
		}
	}

void GSPAligner::SetMx()
	{
	const uint QL = m_Q->GetSeqLength();
	const uint TL = m_T->GetSeqLength();
	Alloc(QL, TL);
	for (uint QPos = 0; QPos < QL; ++QPos)
		{
		for (uint TPos = 0; TPos < TL; ++TPos)
			{
			double Score = GetScore_Pos(QPos, TPos);
			m_Mx[QPos][TPos] = (float) Score;
			}
		}
	}

void GSPAligner::Align(const GSProf &Q, const GSProf &T)
	{
	m_Q = &Q;
	m_T = &T;

	const uint QNF = m_Q->GetFeatureCount();
	const uint TNF = m_T->GetFeatureCount();
	asserta(QNF == TNF);
	m_FeatureCount = QNF;

	const uint QL = Q.GetSeqLength();
	const uint TL = T.GetSeqLength();

	SetMx();
	MxToTsv(g_ftsv);

	string Path;

	XDPMem Mem;
	ViterbiFastMem(Mem, m_Mx, QL, TL, Path);

	const byte *SeqQ = (const byte *) Q.m_Seq.c_str();
	const byte *SeqT = (const byte *) T.m_Seq.c_str();
	LogAlnPretty(SeqQ, SeqT, Path.c_str(), false);
//	Log("Path=%s\n", Path.c_str());
	}
