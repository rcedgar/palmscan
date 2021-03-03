#include "myutils.h"
#include "pssm.h"
#include "sfasta.h"
#include "pssmsearch.h"
#include <set>

unsigned SLCSet(const set<string> &Seqs, map<string, unsigned> &SeqToClusterIndex);

static vector<set<string> > g_UniqueSeqs;
static vector<map<string, unsigned> > g_SeqToClusterIndex;

static void OnHit2(PSSMHitPairNt &Hit)
	{
	const unsigned L = Hit.AASegLen;
	if (L < opt_min_cdr3 || L > opt_max_cdr3)
		return;

	string Seq;
	for (unsigned i = 0; i < L; ++i)
		{
		char c = Hit.AASeg[i];
		if (g_CharToLetterAmino[c] >= 20)
			return;
		Seq.push_back(c);
		}

	asserta(L < SIZE(g_SeqToClusterIndex));
	map<string, unsigned>::const_iterator p = g_SeqToClusterIndex[L].find(Seq);
	asserta(p != g_SeqToClusterIndex[L].end());
	unsigned ClusterIndex = p->second;

	fprintf(g_fOut, "%u", L);
	fprintf(g_fOut, "\t%u", ClusterIndex);
	fprintf(g_fOut, "\t%s", Seq.c_str());
	fprintf(g_fOut, "\t%s", Hit.Label->c_str());
	fprintf(g_fOut, "\n");
	}

static void OnHit(PSSMHitPairNt &Hit)
	{
	const unsigned L = Hit.AASegLen;
	if (L < opt_min_cdr3 || L > opt_max_cdr3)
		return;

	string Seq;
	for (unsigned i = 0; i < L; ++i)
		{
		char c = Hit.AASeg[i];
		if (g_CharToLetterAmino[c] >= 20)
			return;
		Seq.push_back(c);
		}
	g_UniqueSeqs[Hit.AASegLen].insert(Seq);
	}

void SLC()
	{
	asserta(optset_psm1 && optset_psm2);
	const string InputFileName = string(opt_slc);
	const string OutputFileName = string(opt_output);
	asserta(InputFileName != "");
	asserta(OutputFileName != "");
	asserta(opt_min_cdr3 <= opt_max_cdr3);

	g_UniqueSeqs.resize(opt_max_cdr3+1);
	g_SeqToClusterIndex.resize(opt_max_cdr3+1);


	g_fOut = CreateStdioFile(OutputFileName);

	PSSM P1;
	PSSM P2;
	P1.FromFile(opt_psm1);
	P2.FromFile(opt_psm2);

	SFasta SF;
	SF.Open(InputFileName);

	DoSearchPairNt(P1, (float) opt_minscore1, P2, (float) opt_minscore2, SF, OnHit);

	Log("\n");
	Log("Len  Unique seqs\n");
	Log("---  -----------\n");
	for (unsigned L = opt_min_cdr3; L <= opt_max_cdr3; ++L)
		{
		const set<string> &Seqs = g_UniqueSeqs[L];
		unsigned N = SIZE(Seqs);
		Log("%3u  %11u\n", L, N);
		}

	Log("\n");
	unsigned LengthCount = opt_max_cdr3 - opt_min_cdr3 + 1;
	for (unsigned L = opt_min_cdr3; L <= opt_max_cdr3; ++L)
		{
		ProgressStep(L - opt_min_cdr3, LengthCount, "Clustering");
		const set<string> &Seqs = g_UniqueSeqs[L];
		unsigned N = SIZE(Seqs);

		SLCSet(Seqs, g_SeqToClusterIndex[L]);
		asserta(SIZE(g_SeqToClusterIndex[L]) == N);
#if 0
		for (map<string, unsigned>::const_iterator p = SeqToClusterIndex.begin();
		  p != SeqToClusterIndex.end(); ++p)
			{
			const string &Seq = p->first;
			unsigned ClusterIndex = p->second;
			Log("%3u  %10u  %s\n", L, ClusterIndex, Seq.c_str());
			}
#endif
		}

	SF.Rewind();
	ProgressStep(0, 1002, "Writing output");
	for (;;)
		{
		const char *NtSeq = SF.GetNextSeq();
		if (NtSeq == 0)
			break;
		ProgressStep(SF.GetPctDoneX10(), 1002, "Writing output");

		const string Label = SF.GetLabel();
		const unsigned Lnt = SF.GetSeqLength();
		
		PSSMHitPairNt Hit;
		GetTopHitNrPair(P1, (float) opt_minscore1, P2, (float) opt_minscore2, Label, NtSeq, Lnt, Hit);
		OnHit2(Hit);
		}
	ProgressStep(1001, 1002, "Writing output");

	CloseStdioFile(g_fOut);
	}
