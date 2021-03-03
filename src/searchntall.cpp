#include "myutils.h"
#include "sfasta.h"
#include "pssm.h"
#include "pssmsearch.h"

void RevCompInPlace(byte *Seq, unsigned L);

static void SearchNtAll1(FILE *fOut, const PSSM &P, float MinScore, 
  const string &Label, const char *Seq, unsigned L, bool Strand)
	{
	const unsigned M = P.GetColCount();
	if (M > L)
		return;
	float BestScore = -9e9f;
	unsigned BestPos = UINT_MAX;
	for (unsigned Pos = 0; Pos <= L - M; ++Pos)
		{
		float Score = P.GetScore(Seq + Pos);
		if (Score > MinScore)
			{
			fprintf(fOut, "%s", Label.c_str());
			fprintf(fOut, "\t%c", pom(Strand));
			fprintf(fOut, "\t%u", Pos + 1);
			fprintf(fOut, "\t%.3f", Score);
			fprintf(fOut, "\t%*.*s", M, M, Seq + Pos);
			fprintf(fOut, "\n");
			}
		}
	}

void SearchNtAll()
	{
	asserta(optset_psm);
	const string InputFileName = string(opt_search_nt_all);
	const string OutputFileName = string(opt_output);
	asserta(InputFileName != "");
	asserta(OutputFileName != "");

	FILE *fOut = CreateStdioFile(OutputFileName);

	PSSM P;
	P.FromFile(opt_psm);
	unsigned M = P.GetColCount();

	SFasta SF;
	SF.Open(InputFileName);

	for (;;)
		{
		const char *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;

		const string Label = SF.GetLabel();
		const unsigned L = SF.GetSeqLength();
		SearchNtAll1(fOut, P, (float) opt_minscore, Label, Seq, L, true);

		RevCompInPlace((byte *) Seq, L);

		SearchNtAll1(fOut, P, (float) opt_minscore, Label, Seq, L, false);
		}

	CloseStdioFile(fOut);
	}
