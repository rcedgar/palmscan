#include "myutils.h"
#include "tshitmgr.h"
#include "sort.h"
#include "abcxyz.h"
#include "omplock.h"
#include "outputfiles.h"

double GetConfidenceScore(double MotifRMSD);

void TSHitMgr::SetQuery(PDBChain &Query)
	{
	m_Hits.clear();
	m_Query = &Query;
	}

void TSHitMgr::Add(TSHit &TH)
	{
	m_Hits.push_back(TH);
	}

void TSHitMgr::SetTopHit()
	{
	double BestMotifRMSD2 = DBL_MAX;
	TSHit *BestHit = 0;
	for (uint i = 0; i < SIZE(m_Hits); ++i)
		{
		TSHit &Hit = m_Hits[i];
		if (Hit.m_MotifRMSD2 < BestMotifRMSD2)
			{
			BestMotifRMSD2 = Hit.m_MotifRMSD2;
			BestHit = &Hit;
			}
		}
	m_TopHit = BestHit;
	}

void TSHitMgr::WriteReport(FILE *f) const
	{
	if (f == 0)
		return;

	asserta(m_Query != 0);

	const string &QLabel = m_Query->m_Label;
	uint QSeqLength = m_Query->GetSeqLength();

	fprintf(f, "\n");
	fprintf(f, "_______________________________________________\n");
	fprintf(f, "Query >%s (%u aa)\n", QLabel.c_str(), QSeqLength);

	vector<double> Ranks;
	const uint N = SIZE(m_Hits);
	if (N == 0)
		return;
	for (uint i = 0; i < N; ++i)
		{
		const TSHit &Hit = m_Hits[i];
		Ranks.push_back(Hit.m_MotifRMSD2);
		}

	vector<uint> Order(N);
	QuickSortOrder(Ranks.data(), N, Order.data());

	fprintf(f, "%5.5s", "Score");
	fprintf(f, " %6.6s", "RMSDm");
	fprintf(f, " %5.5s", "PosA");
	fprintf(f, " %5.5s", "PosB");
	fprintf(f, " %5.5s", "PosC");
	fprintf(f, "  %12.12s", "MotifA");
	fprintf(f, "  %14.14s", "MotifB");
	fprintf(f, "  %8.8s", "MotifC");
	fprintf(f, "  %s", "Ref");
	fprintf(f, "\n");

	uint n = min(N, 4u);
	double TopScore = 0;
	for (uint i = 0; i < n; ++i)
		{
		uint k = Order[i];

		const TSHit &TH = m_Hits[k];

		uint QPosA = TH.m_QPosA;
		uint QPosB = TH.m_QPosB;
		uint QPosC = TH.m_QPosC;
		double MotifRMSD2 = TH.m_MotifRMSD2;
		double MotifRMSD = sqrt(MotifRMSD2);
		double ConfScore = GetConfidenceScore(MotifRMSD);
		if (i == 0)
			TopScore = ConfScore;

		const string &RefLabel = TH.m_Ref->m_Label;

		string ASeq, BSeq, CSeq;
		TH.m_Query->GetSubSeq(QPosA, AL, true, ASeq);
		TH.m_Query->GetSubSeq(QPosB, BL, true, BSeq);
		TH.m_Query->GetSubSeq(QPosC, CL, true, CSeq);

		fprintf(f, "%5.2f", ConfScore);
		fprintf(f, " %6.2f", MotifRMSD);
		fprintf(f, " %5u", QPosA + 1);
		fprintf(f, " %5u", QPosB + 1);
		fprintf(f, " %5u", QPosC + 1);
		fprintf(f, "  %s", ASeq.c_str());
		fprintf(f, "  %s", BSeq.c_str());
		fprintf(f, "  %s", CSeq.c_str());
		fprintf(f, "  %s", RefLabel.c_str());
		fprintf(f, "\n");
		}

	double SketchPct = m_TopHit->m_SketchPct;
	m_TopHit->WriteSketch(f);
	m_TopHit->WriteAln(f);

	double FinalScore = m_TopHit->m_Score;

	fprintf(f, "\n");
	fprintf(f, "Score %.3f (%.2f, %.1f%%)", FinalScore, TopScore, SketchPct);
	if (TopScore >= 0.75)
		fprintf(f, " high confidence");
	else if (TopScore >= 0.50)
		fprintf(f, " low confidence");
	else
		fprintf(f, " likely false positive");
	fprintf(f, "\n");
	}

void TSHitMgr::WriteOutput()
	{
	LockOutput();
	SetTopHit();
	if (m_TopHit == 0)
		{
		UnlockOutput();
		return;
		}
	m_TopHit->SetSketch();
	m_TopHit->WriteTsv(g_ftsv);
	m_TopHit->WritePalmprintFasta(g_fppfa);
	WritePPC(g_fppc);
	WriteReport(g_freport_3d);
	UnlockOutput();
	}

void TSHitMgr::GetPPC(PDBChain &PPC) const
	{
	asserta(m_TopHit != 0);
	asserta(m_Query != 0);

	const PDBChain &Q = *m_Query;

	uint QPosA = m_TopHit->m_QPosA;
	uint QPosB = m_TopHit->m_QPosB;
	uint QPosC = m_TopHit->m_QPosC;

	asserta(QPosA < QPosB&& QPosB < QPosC);
	uint PPL = QPosC + CL - QPosA;

	vector<vector<double> > MotifCoords(3);
	Q.GetPt(QPosA, MotifCoords[A]);
	Q.GetPt(QPosB, MotifCoords[B]);
	Q.GetPt(QPosC, MotifCoords[C]);

	vector<double> t;
	vector<vector<double> > R;
	GetTriForm(MotifCoords, t, R);

	PPC.Clear();
	PPC.m_Label = Q.m_Label;
	PPC.m_MotifPosVec.clear();
	PPC.m_MotifPosVec.push_back(0);
	PPC.m_MotifPosVec.push_back(QPosB - QPosA);
	PPC.m_MotifPosVec.push_back(QPosC - QPosA);

	const uint N = SIZE(Q.m_Xs);
	asserta(SIZE(Q.m_Seq) == N);
	asserta(SIZE(Q.m_Ys) == N);
	asserta(SIZE(Q.m_Zs) == N);
	vector<double> Pt(3);
	vector<double> XPt(3);
	for (uint i = 0; i < PPL; ++i)
		{
		char c = Q.m_Seq[QPosA + i];
		PPC.m_Seq += c;

		Q.GetPt(QPosA + i, Pt);
		XFormPt(Pt, t, R, XPt);
		PPC.m_Xs.push_back(XPt[X]);
		PPC.m_Ys.push_back(XPt[Y]);
		PPC.m_Zs.push_back(XPt[Z]);
		}
	}

void TSHitMgr::WritePPC(FILE* f) const
	{
	if (f == 0)
		return;
	if (m_TopHit == 0)
		return;

	uint PosA = m_TopHit->m_QPosA;
	uint PosB = m_TopHit->m_QPosB;
	uint PosC = m_TopHit->m_QPosC;

	PDBChain PPC;
	m_Query->SetMotifPosVec(PosA, PosB, PosC);
	m_Query->GetPPC(PPC);
	asserta(PPC.CheckPPCMotifCoords());
	PPC.ToCal(f);
	}
