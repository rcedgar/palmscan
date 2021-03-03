#include "myutils.h"
#include "pssm.h"
#include "sfasta.h"
#include "pssmsearch.h"
#include <set>

unsigned SLCSet(const set<string> &Seqs, map<string, unsigned> &SeqToClusterIndex);

static vector<set<string> > g_UniqueSeqs;
static vector<map<string, unsigned> > g_SeqToClusterIndex;

void SLC2()
	{
	const string InputFileName = string(opt_slc2);
	const string OutputFileName = string(opt_output);
	asserta(InputFileName != "");
	asserta(OutputFileName != "");
	asserta(opt_min_cdr3 <= opt_max_cdr3);

	FILE *fIn = OpenStdioFile(InputFileName);
	FILE *fOut = CreateStdioFile(OutputFileName);

	g_UniqueSeqs.resize(opt_max_cdr3+1);
	g_SeqToClusterIndex.resize(opt_max_cdr3+1);

	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		const string &Label = Fields[0];
		const string &CDR3 = Fields[1];
		const unsigned L = SIZE(CDR3);
		if (L < opt_min_cdr3 || L > opt_max_cdr3)
			continue;

		g_UniqueSeqs[L].insert(CDR3);
		}

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
	unsigned Total = 0;
	for (unsigned L = opt_min_cdr3; L <= opt_max_cdr3; ++L)
		{
		ProgressStep(L - opt_min_cdr3, LengthCount, "Clustering");
		const set<string> &Seqs = g_UniqueSeqs[L];
		unsigned N = SIZE(Seqs);

		unsigned ClusterCount = SLCSet(Seqs, g_SeqToClusterIndex[L]);
		asserta(SIZE(g_SeqToClusterIndex[L]) == N);
		for (map<string, unsigned>::iterator p = g_SeqToClusterIndex[L].begin();
		  p != g_SeqToClusterIndex[L].end(); ++p)
			 p->second += Total;
		Total += ClusterCount;

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

	SetStdioFilePos(fIn, 0);
	while (ReadLineStdioFile(fIn, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);
		const string &Label = Fields[0];
		const string &CDR3 = Fields[1];
		const unsigned L = SIZE(CDR3);
		if (L < opt_min_cdr3 || L > opt_max_cdr3)
			continue;

		map<string, unsigned>::const_iterator p = g_SeqToClusterIndex[L].find(CDR3);
		asserta(p != g_SeqToClusterIndex[L].end());
		unsigned ClusterIndex = p->second;
		fprintf(fOut, "%u", ClusterIndex);
		fprintf(fOut, "\t%s", Label.c_str());
		fprintf(fOut, "\t%s", CDR3.c_str());
		fputc('\n', fOut);
		}
	ProgressStep(1001, 1002, "Writing output");

	CloseStdioFile(fIn);
	CloseStdioFile(fOut);
	}
