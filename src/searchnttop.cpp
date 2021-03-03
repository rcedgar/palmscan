#include "myutils.h"
#include "sfasta.h"
#include "pssm.h"
#include "pssmsearch.h"

static unsigned Xlat(const char *nt, unsigned Lnt, char *aa)
	{
	Lnt -= Lnt%3;
	unsigned Laa = 0;
	for (unsigned Pos = 0; Pos < Lnt; Pos += 3)
		{
		byte c1 = nt[Pos];
		byte c2 = nt[Pos+1];
		byte c3 = nt[Pos+2];

		byte x1 = g_CharToLetterNucleo[c1];
		byte x2 = g_CharToLetterNucleo[c2];
		byte x3 = g_CharToLetterNucleo[c3];

		unsigned Word = 16*x1 + 4*x2 + x3;
		byte AminoChar = 'X';
		if (Word < 64)
			AminoChar = g_CodonWordToAminoChar[Word];
		aa[Laa++] = AminoChar;
		}
	//Log("%*.*s\n", Laa, Laa, aa);
	return Laa;
	}

void SearchNtTop()
	{
	asserta(optset_psm);
	const string InputFileName = string(opt_search_nt_top);
	const string OutputFileName = string(opt_output);
	asserta(InputFileName != "");
	asserta(OutputFileName != "");

	FILE *fOut = CreateStdioFile(OutputFileName);

	PSSM P;
	P.FromFile(opt_psm);
	unsigned M = P.GetColCount();

	SFasta SF;
	SF.Open(InputFileName);

	char *AASeq = 0;
	unsigned MaxLaa = 0;
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;

		const string Label = SF.GetLabel();
		const unsigned Lnt = SF.GetSeqLength();

		if ((Lnt+2)/3 > MaxLaa)
			{
			if (AASeq != 0)
				myfree(AASeq);
			MaxLaa = Lnt;
			AASeq = myalloc(char, MaxLaa);
			}

		unsigned Laa = Xlat(Seq, Lnt, AASeq);

		PSSMHitAA Hit;
		GetTopHitAA(P, (float) opt_minscore, Label, (const char *) AASeq, Laa, Hit);
		unsigned AAPos = Hit.AAPos;
		if (AAPos == UINT_MAX)
			continue;

		fprintf(fOut, "%s", Label.c_str());
		fprintf(fOut, "\t%u", AAPos + 1);
		fprintf(fOut, "\t%.3f", Hit.Score);
		fprintf(fOut, "\t%*.*s", M, M, AASeq + AAPos);
		fprintf(fOut, "\t%u", 3*AAPos + 1);
		fprintf(fOut, "\t%*.*s", 3*M, 3*M, Seq + 3*AAPos);
		fprintf(fOut, "\n");
		}

	CloseStdioFile(fOut);
	}
