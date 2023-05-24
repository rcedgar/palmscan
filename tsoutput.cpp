#include "myutils.h"
#include "trisearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "omplock.h"

double GetConfidenceScore(double MotifRMSD);

void TriSearcher::WriteOutput()
	{
	LockOutput();
	WriteTsv();
	WriteReport();
	UnlockOutput();
	}

void TriSearcher::WriteTsv()
	{
	FILE *f = g_ftsv_tri;
	if (f == 0)
		return;

	TSHit TH;
	bool Ok = GetTopHit(TH);
	if (Ok)
		TH.WriteTsv(f);
	}

void TriSearcher::WriteReport()
	{
	FILE *f = g_freport_tri;
	if (f == 0)
		return;
	const uint TriHitCount = SIZE(m_HitOrder);
	if (TriHitCount == 0)
		return;


	const string &QueryLabel = m_Query->m_Label;
	const string &RefLabel = m_Ref->m_Label;

	const char *QLabel = QueryLabel.c_str();
	const char *RLabel = RefLabel.c_str();
	uint QSeqLength = SIZE(m_Query->m_Seq);
	uint RSeqLength = SIZE(m_Ref->m_Seq);

	fprintf(f, "\n");
	fprintf(f, "_______________________________________________\n");
	fprintf(f, "Query >%s (%u aa)\n", QLabel, QSeqLength);
	fprintf(f, "  Ref >%s (%u aa)", RLabel, RSeqLength);
	fprintf(f, ", %u tris\n", TriHitCount);

	fprintf(f, "%5.5s", "Score");
	fprintf(f, "  %6.6s", "RMSDm");
	fprintf(f, "  %6.6s", "RMSDt");
	fprintf(f, " %5.5s", "PosA");
	fprintf(f, " %5.5s", "PosB");
	fprintf(f, " %5.5s", "PosC");
	fprintf(f, "  %12.12s", "MotifA");
	fprintf(f, "  %14.14s", "MotifB");
	fprintf(f, "  %8.8s", "MotifC");
	fprintf(f, "\n");

	uint n = min(TriHitCount, 4u);
	for (uint i = 0; i < n; ++i)
		{
		uint k = m_HitOrder[i];

		uint PosA = m_PosAs[k];
		uint PosB = m_PosBs[k];
		uint PosC = m_PosCs[k];
		double TriRMSD2 = m_TriRMSD2s[k];
		double MotifRMSD2 = m_MotifRMSD2s[k];

		double MotifRMSD = sqrt(MotifRMSD2);
		double TriRMSD = sqrt(TriRMSD2);
		double ConfScore = GetConfidenceScore(MotifRMSD);

		string ASeq, BSeq, CSeq;
		m_Query->GetSubSeq(PosA, AL, true, ASeq);
		m_Query->GetSubSeq(PosB, BL, true, BSeq);
		m_Query->GetSubSeq(PosC, CL, true, CSeq);

		fprintf(f, "%5.2f", ConfScore);
		fprintf(f, "  %6.2f", MotifRMSD);
		fprintf(f, "  %6.2f", TriRMSD);
		fprintf(f, " %5u", PosA + 1);
		fprintf(f, " %5u", PosB + 1);
		fprintf(f, " %5u", PosC + 1);
		fprintf(f, "  %s", ASeq.c_str());
		fprintf(f, "  %s", BSeq.c_str());
		fprintf(f, "  %s", CSeq.c_str());
		fprintf(f, "\n");
		}

	WriteAln(f);

	TSHit TopHit;
	bool Ok = GetTopHit(TopHit);
	if (Ok)
		{
		TopHit.SetSketch();
		TopHit.WriteSketch(f);
		}
	}

void TriSearcher::WriteAln(FILE *f)
	{
	if (f == 0)
		return;

	TSHit TopHit;
	bool Ok = GetTopHit(TopHit);
	if (!Ok)
		return;
	TopHit.WriteAln(f);
	}
