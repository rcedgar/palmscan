#include "myutils.h"
#include "pssm.h"
#include "sfasta.h"
#include "pssmsearch.h"
#include "sort.h"
#include <set>

// Input is output from assignv: Query TargetV PctId
// CDR3 sequence is in label.
// readid=HV2697C05FUMVM;spec=113;stage=Acute;repl=A;cdr3=AKFGYSSTWFTQGLDV;        IGHV1-69     89.5

unsigned SLCSet(const set<string> &Seqs, map<string, unsigned> &SeqToClusterIndex);
double MakeClusterStr(unsigned L, unsigned ClusterIndex, const vector<string> &Labels,
  string &s);

static vector<set<string> > g_UniqueSeqs;
static vector<map<string, unsigned> > g_SeqToClusterIndex;
static map<pair<unsigned, unsigned>, vector<string> > g_LengthAndClusterIndexToLabels;

// readid=HV2697C05FUMVM;spec=113;stage=Acute;repl=A;cdr3=AKFGYSSTWFTQGLDV;
void GetCDR3_aa(const string &Label, string &CDR3_aa)
	{
	CDR3_aa.clear();
	const char *s = Label.c_str();
	const char *s1 = strstr(s, "cdr3=");
	if (s1 == 0)
		Die("cdr3= not found in >%s", s);
	s1 += 5;
	const char *s2 = strchr(s1, ';');
	if (s2 == 0)
		Die("Missing ; after cdr3= in >%s", s);

	for (const char *p = s1; p != s2; ++p)
		{
		char c = *p;
		if (c == '*')
			{
			CDR3_aa = "(stop)";
			return;
			}
		CDR3_aa.push_back(c);
		}
	}

static void Pass2(const string &Label)
	{
	string CDR3_aa;
	GetCDR3_aa(Label, CDR3_aa);
	asserta(!CDR3_aa.empty());
	if (CDR3_aa[0] == '(')
		return;

	unsigned L = SIZE(CDR3_aa);
	if (L < opt_min_cdr3 || L > opt_max_cdr3)
		return;

	asserta(L < SIZE(g_SeqToClusterIndex));
	map<string, unsigned>::const_iterator p = g_SeqToClusterIndex[L].find(CDR3_aa);
	asserta(p != g_SeqToClusterIndex[L].end());
	unsigned ClusterIndex = p->second;
	pair<unsigned, unsigned> LengthAndClusterIndex(L, ClusterIndex);
	map<pair<unsigned, unsigned>, vector<string> >::const_iterator q =
	  g_LengthAndClusterIndexToLabels.find(LengthAndClusterIndex);
	if (q == g_LengthAndClusterIndexToLabels.end())
		{
		vector<string> Empty;
		g_LengthAndClusterIndexToLabels[LengthAndClusterIndex] = Empty;
		}
	g_LengthAndClusterIndexToLabels[LengthAndClusterIndex].push_back(Label);

	//fprintf(g_fOut, "%u", L);
	//fprintf(g_fOut, "\t%u", ClusterIndex);
	//fprintf(g_fOut, "\t%s", CDR3_aa.c_str());
	//fprintf(g_fOut, "\t%s", Label.c_str());
	//fprintf(g_fOut, "\n");
	}

static void Pass1(const string &Label)
	{
	string CDR3_aa;
	GetCDR3_aa(Label, CDR3_aa);
	asserta(!CDR3_aa.empty());
	if (CDR3_aa[0] == '(')
		return;

	unsigned L = SIZE(CDR3_aa);
	if (L < opt_min_cdr3 || L > opt_max_cdr3)
		return;

	g_UniqueSeqs[L].insert(CDR3_aa);
	}

void SLCU()
	{
	const string InputFileName = string(opt_slcu);
	const string OutputFileName = string(opt_output);
	const string ReportFileName = string(opt_report);
	asserta(InputFileName != "");
	asserta(OutputFileName != "");
	asserta(opt_min_cdr3 <= opt_max_cdr3);

	g_UniqueSeqs.resize(opt_max_cdr3+1);
	g_SeqToClusterIndex.resize(opt_max_cdr3+1);

	FILE *fIn = OpenStdioFile(InputFileName);
	g_fOut = CreateStdioFile(OutputFileName);
	if (ReportFileName != "")
		g_fRep = CreateStdioFile(ReportFileName);

	string Line;
	vector<string> Fields;
	ProgressFileInit(fIn, "SLCU pass 1 %s", InputFileName.c_str());
	while (ReadLineStdioFile(fIn, Line))
		{
		ProgressFileStep();
		Pass1(Line);
		}
	ProgressFileDone();

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

	SetStdioFilePos64(fIn, 0);
	ProgressFileInit(fIn, "SLCU pass 2 %s", InputFileName.c_str());
	while (ReadLineStdioFile(fIn, Line))
		{
		ProgressFileStep();
		Split(Line, Fields, '\t');
		Pass2(Line);
		}
	ProgressFileDone();

	vector<float> Scores;
	vector<string> Recs;
	unsigned ClusterIndex2 = 0;
	unsigned ClusterCount = SIZE(g_LengthAndClusterIndexToLabels);
	for (map<pair<unsigned, unsigned>, vector<string> >::const_iterator p =
	  g_LengthAndClusterIndexToLabels.begin(); p != g_LengthAndClusterIndexToLabels.end();
	  ++p)
		{
		ProgressStep(ClusterIndex2, ClusterCount, "Writing output");
		pair<unsigned, unsigned> LengthAndClusterIndex = p->first;
		unsigned L = LengthAndClusterIndex.first;
		const vector<string> &Labels = p->second;
		string Rec;
		double Score = MakeClusterStr(L, ClusterIndex2, Labels, Rec);
		Scores.push_back((float) Score);
		Recs.push_back(Rec);
		const unsigned N = SIZE(Labels);
		asserta(N > 0);
		for (unsigned i = 0; i < N; ++i)
			fprintf(g_fOut, "%u\t%s\n", ClusterIndex2, Labels[i].c_str());
		++ClusterIndex2;
		}

	if (g_fRep != 0)
		{
		vector<unsigned> Order;
		asserta(SIZE(Scores) == ClusterCount);
		SortDescending(Scores, Order);
		asserta(SIZE(Order) == ClusterCount);
		for (unsigned i = 0; i < ClusterCount; ++i)
			{
			unsigned k = Order[i];
			float Score = Scores[k];
			const string &Rec = Recs[k];
			fprintf(g_fRep, "%.1f\trank=%u\t", Score, i+1);
			fputs(Rec.c_str(), g_fRep);
			fputc('\n', g_fRep);
			}
		}

	CloseStdioFile(fIn);
	CloseStdioFile(g_fOut);
	CloseStdioFile(g_fRep);
	}
