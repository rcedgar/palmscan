#include "myutils.h"
#include "trisearcher.h"
#include "pdbchain.h"
#include "sort.h"
#include "abcxyz.h"
#include "searchparams.h"

const uint MotifLVec[3] = { 12, 14, 8 };

void TriSearcher::LogMe(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "TriSearcher() Radius %.1f\n", Radius);
	fprintf(f, "NAB %u..%u", NABmin, NABmax);
	fprintf(f, ", NBC %u..%u", NBCmin, NBCmax);
	fprintf(f, ", NAC %u..%u\n", NACmin, NACmax);

	const uint N = SIZE(m_TriRMSD2s);
	if (N > 0)
		{
		Log("%u hits\n", N);
		Log(" PosA   PosB   PosC  TriD2    MotifD2\n");
		//   12345  12345  12345  12345  123456789
		if (m_Query != 0 && m_Query->m_MotifPosVec.size() == 3)
			{
			string QA, QB, QC;
			uint QPosA = m_Query->m_MotifPosVec[0];
			uint QPosB = m_Query->m_MotifPosVec[1];
			uint QPosC = m_Query->m_MotifPosVec[2];
			m_Query->GetMotifSeq(0, QA);
			m_Query->GetMotifSeq(1, QB);
			m_Query->GetMotifSeq(2, QC);
			Log("%5u", QPosA);
			Log("  %5u", QPosB);
			Log("  %5u", QPosC);
			Log("  %5.5s", "");
			Log("  %9.9s", "");
			Log("  %s", QA.c_str());
			Log("  %s", QB.c_str());
			Log("  %s", QC.c_str());
			Log("  << Query\n");
			}

		for (uint k = 0; k < N; ++k)
			{
			uint i = m_HitOrder[k];
			uint PosA = m_PosAs[i];
			uint PosB = m_PosBs[i];
			uint PosC = m_PosCs[i];
			double TriRMSD2 = m_TriRMSD2s[i];
			double MotifRMSD2 = m_MotifRMSD2s[i];

			string A, B, C;
			m_Query->GetSubSeq(PosA, AL, false, A);
			m_Query->GetSubSeq(PosB, BL, false, B);
			m_Query->GetSubSeq(PosC, CL, false, C);

			Log("%5u", PosA);
			Log("  %5u", PosB);
			Log("  %5u", PosC);
			Log("  %5.2f", TriRMSD2);
			Log("  %9.1f", MotifRMSD2);
			Log("  %s", A.c_str());
			Log("  %s", B.c_str());
			Log("  %s", C.c_str());
			if (k == 0)
				{
				double TM = GetMotifsTM();
				Log("  TM=%.3g", TM);
				}
			Log("\n");
			}
		}
	}

void TriSearcher::SetHitOrder()
	{
	const uint N = SIZE(m_PosAs);
	if (N == 0)
		{
		m_HitOrder.clear();
		return;
		}
	m_HitOrder.resize(N);
	QuickSortOrder(m_MotifRMSD2s.data(), N, m_HitOrder.data());
	}

void TriSearcher::SetTriForm()
	{
	asserta(!m_HitOrder.empty());

	uint k = m_HitOrder[0];
	uint QueryPosA = m_PosAs[k];
	uint QueryPosB = m_PosBs[k];
	uint QueryPosC = m_PosCs[k];

	vector<vector<double> > QueryMotifCoords(3);
	m_Query->GetPt(QueryPosA, QueryMotifCoords[A]);
	m_Query->GetPt(QueryPosB, QueryMotifCoords[B]);
	m_Query->GetPt(QueryPosC, QueryMotifCoords[C]);

	GetTriForm(QueryMotifCoords, m_TriForm_t, m_TriForm_R);
	}

double TriSearcher::GetRMSDMotifs(uint QueryPosA, uint QueryPosB,
  uint QueryPosC) const
	{
	vector<vector<double> > RefMotifCoords;
	m_Ref->GetMotifCoords(RefMotifCoords);

	asserta(m_Ref != 0 && m_Ref->m_MotifPosVec.size() == 3);
	uint RefPosA = m_Ref->m_MotifPosVec[A];
	uint RefPosB = m_Ref->m_MotifPosVec[B];
	uint RefPosC = m_Ref->m_MotifPosVec[C];

	vector<vector<double> > QueryMotifCoords(3);
	m_Query->GetPt(QueryPosA, QueryMotifCoords[A]);
	m_Query->GetPt(QueryPosB, QueryMotifCoords[B]);
	m_Query->GetPt(QueryPosC, QueryMotifCoords[C]);

	vector<double> t;
	vector<vector<double> > R;
	GetTriForm(QueryMotifCoords, t, R);

	double RMSD2A = GetRMSD2Segment(QueryPosA, RefPosA, AL, t, R);
	double RMSD2B = GetRMSD2Segment(QueryPosB, RefPosB, BL, t, R);
	double RMSD2C = GetRMSD2Segment(QueryPosC, RefPosC, CL, t, R);
	double RMSD2 = RMSD2A + RMSD2B + RMSD2C;
	double Mean = (RMSD2A + RMSD2B + RMSD2C)/(AL + BL + CL);
	double RMSD = sqrt(Mean);
	return RMSD;
	}

double TriSearcher::GetRMSD2Segment(uint QueryStartPos, uint RefStartPos,
  uint n, const vector<double> &t, const vector<vector<double> > &R) const
	{
	double Sum = 0;

	vector<double> RefPt(3);
	vector<double> QueryPt(3);
	vector<double> QueryPtX(3);
	for (uint i = 0; i < n; ++i)
		{
		m_Ref->GetPt(RefStartPos + i, RefPt);
		m_Query->GetPt(QueryStartPos + i, QueryPt);

		XFormPt(QueryPt, t, R, QueryPtX);

		double d2 = GetDist2(RefPt, QueryPtX);
		Sum += d2;
		}

	return Sum;
	}

bool TriSearcher::GetTopHit(uint &QueryPosA, uint &QueryPosB, uint &QueryPosC) const
	{
	if (m_HitOrder.empty())
		{
		QueryPosA = UINT_MAX;
		QueryPosB = UINT_MAX;
		QueryPosC = UINT_MAX;
		return false;
		}

	uint k = m_HitOrder[0];
	QueryPosA = m_PosAs[k];
	QueryPosB = m_PosBs[k];
	QueryPosC = m_PosCs[k];
	return true;
	}

void TriSearcher::GetHit(uint Index, TSHit &TH) const
	{
	asserta(Index < SIZE(m_HitOrder));
	uint k = m_HitOrder[Index];

	TH.m_Query = m_Query;
	TH.m_Ref = m_Ref;

	TH.m_MotifRMSD2 = m_MotifRMSD2s[k];
	TH.m_TriRMSD2 = m_TriRMSD2s[k];

	TH.m_QPosA = m_PosAs[k];
	TH.m_QPosB = m_PosBs[k];
	TH.m_QPosC = m_PosCs[k];

	asserta(m_Ref != 0 && m_Ref->m_MotifPosVec.size() == 3);
	}

bool TriSearcher::GetTopHit(TSHit &TH) const
	{
	if (m_HitOrder.empty())
		{
		TH.Clear();
		return false;
		}
	GetHit(0, TH);
	return true;
	}
