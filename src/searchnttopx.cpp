#include "myutils.h"
#include "sfasta.h"
#include "xlat.h"
#include "pssm.h"
#include "pssmsearch.h"

#define TRACE	1

extern int g_Frame;

void GetTopHitX(const PSSM &P, float MinScore,
  const string &Label, const char *NtSeq, unsigned Lnt,
  PSSMHitNt &Hit)
	{
#if	TRACE
	Log("\n");
	Log("GetTopHixX %umt >%s\n", Lnt, Label.c_str());
	Log("%*.*s\n", Lnt, Lnt, NtSeq);
#endif

	float BestScore = -9e9f;
	unsigned BestAAPos = UINT_MAX;
	int BestFrame = 0;

	AllocXlat(Lnt);

	for (int Frame = -3; Frame <= +3; ++Frame)
		{
		if (Frame == 0)
			continue;
		if (g_Frame != 0 && Frame != g_Frame)
			continue;

		byte *AASeq = g_AASeqs[Frame+3];
		unsigned Laa = XlatX((const byte *) NtSeq, Lnt, Frame, AASeq);
#if	TRACE
		Log("%+2d %*.*s\n", Frame, Laa, Laa, AASeq);
#endif

		PSSMHitAA Hit2;
		GetTopHitAA(P, MinScore, Label, (const char *) AASeq, Laa, Hit2);
		if (Hit2.AAPos == UINT_MAX)
			{
#if	TRACE
			Log("Hit2.AAPos=UINT_MAX\n");
#endif
			continue;
			}
#if	TRACE
		Log("Score %.1f, AAPos %u\n", Hit2.Score, Hit2.AAPos);
#endif

#if	DEBUG
		{
		unsigned NtPos = AAPosToNtPos(Hit2.AAPos, Lnt, Frame);
		unsigned M = P.GetColCount();
		asserta(NtPos + 3*M <= Lnt);

		string aa;
		XlatStr(NtSeq, Lnt, NtPos, M, Frame, aa);
		for (unsigned i = 0; i < M; ++i)
			asserta(AASeq[Hit2.AAPos + i] == aa[i]);
		}
#endif
		if (Hit2.Score > BestScore)
			{
			BestScore = Hit2.Score;
			BestAAPos = Hit2.AAPos;
			BestFrame = Frame;
			}
		}

	unsigned NtPos = UINT_MAX;
	if (BestAAPos != UINT_MAX)
		{
		NtPos = AAPosToNtPos(BestAAPos, Lnt, BestFrame);
		unsigned M = P.GetColCount();
		asserta(NtPos + 3*M <= Lnt);
		}

	Hit.Label = &Label;
	Hit.NtSeq = NtSeq;
	Hit.Lnt = Lnt;
	Hit.AAPos = BestAAPos;
	Hit.Frame = BestFrame;
	Hit.Score = BestScore;
	Hit.NtPos = NtPos;
	}

void SearchNtTopX()
	{
	asserta(optset_psm);
	const string InputFileName = string(opt_search_nt_topx);
	const string OutputFileName = string(opt_output);
	asserta(InputFileName != "");
	asserta(OutputFileName != "");

	FILE *fOut = CreateStdioFile(OutputFileName);

	PSSM P;
	P.FromFile(opt_psm);
	unsigned M = P.GetColCount();

	SFasta SF;
	SF.Open(InputFileName);

	ProgressStep(0, 1002, "Search topx %s", InputFileName.c_str());
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		ProgressStep(SF.GetPctDoneX10(), 1002, "Search topx %s", InputFileName.c_str());

		const string Label = SF.GetLabel();
		const unsigned Lnt = SF.GetSeqLength();

		PSSMHitNt Hit;
		GetTopHitX(P, (float) opt_minscore, Label, Seq, Lnt, Hit);
		if (Hit.AAPos == UINT_MAX)
			continue;

		fprintf(fOut, "%s", Label.c_str());
		fprintf(fOut, "\t%u", Hit.NtPos);
		fprintf(fOut, "\t%+d", Hit.Frame);
		fprintf(fOut, "\t%.3f", Hit.Score);
		fprintf(fOut, "\t%*.*s", M, M, g_AASeqs[Hit.Frame+3] + Hit.AAPos);
		fprintf(fOut, "\n");
		}
	ProgressStep(1001, 1002, "Search topx %s", InputFileName.c_str());

	CloseStdioFile(fOut);
	}
