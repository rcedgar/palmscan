#include "myutils.h"
#include "sfasta.h"
#include "pssm.h"
#include "pssmsearch.h"

#define TRACE	0

void GetTopHitAA(const PSSM &P, float MinScore, 
  const string &Label, const char *AASeq, unsigned L,
  PSSMHitAA &Hit)
	{
	Hit.Label = &Label;
	Hit.AASeq = AASeq;
	Hit.L = L;
	Hit.Score = -9e9f;
	Hit.AAPos = UINT_MAX;
	asserta(!P.m_IsNucleo);
	const unsigned M = P.GetColCount();
	if (M > L)
		return;
	float BestScore = -9e9f;
	unsigned BestPos = UINT_MAX;
	for (unsigned Pos2 = 0; Pos2 <= L - M; ++Pos2)
		{
		float Score2 = P.GetScore(AASeq + Pos2);
		if (Score2 > BestScore)
			{
			BestScore = Score2;
			BestPos = Pos2;
			}
		}
	Hit.Score = BestScore;
	Hit.AAPos = BestPos;
#if	TRACE
	Log("GetTopHitAA Score %.2f, MinScore %.2f, Pos %u\n",
	  BestScore, MinScore, BestPos);
#endif
	}

void SearchAATop()
	{
	asserta(optset_psm);
	const string InputFileName = string(opt_search_aa_top);
	const string OutputFileName = string(opt_output);
	asserta(InputFileName != "");
	asserta(OutputFileName != "");

	FILE *fOut = CreateStdioFile(OutputFileName);

	PSSM P;
	P.FromFile(opt_psm);
	unsigned M = P.GetColCount();

	SFasta SF;
	SF.Open(InputFileName);

	ProgressStep(0, 1002, "Search aa %s", SF.m_FileName.c_str());
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		ProgressStep(SF.GetPctDoneX10(), 1002, "Search aa %s", SF.m_FileName.c_str());

		const string Label = SF.GetLabel();
		const unsigned L = SF.GetSeqLength();

		PSSMHitAA Hit;
		GetTopHitAA(P, (float) opt_minscore, Label, (const char *) Seq, L, Hit);
		unsigned Pos = Hit.AAPos;
		if (Pos == UINT_MAX)
			continue;

		fprintf(fOut, "%s", Label.c_str());
		fprintf(fOut, "\t%u", Pos + 1);
		fprintf(fOut, "\t%.3f", Hit.Score);
		fprintf(fOut, "\t%*.*s", M, M, Seq + Pos);
		fprintf(fOut, "\n");
		}
	ProgressStep(1001, 1002, "Search aa %s", SF.m_FileName.c_str());

	CloseStdioFile(fOut);
	}
