#include "myutils.h"
#include "sort.h"
#include <map>
#include <set>

// @@ BUG >D00004;taxid=9606;

/***
otus.tsv
--------
Accession       Source			U       S       G       F       Taxon
QOQ00365        gb241_aa        U12391  S2525   G1885   f1343   2697049
KX453308        gb241_nt        .       .       .       .       208893
KX453309        gb241_nt        .       .       .       .       208893
QKO24379        gb241_aa        U37700  s10635  g7759   f1715   11320

taxon.tsv
---------
287144  superkingdom:Viruses,clade:Riboviria,kingdom:Orthornavirae,phylum:Negarnaviricota,subphylum
287145  superkingdom:Viruses,clade:Riboviria,kingdom:Orthornavirae,phylum:Negarnaviricota,subphylum
1774518 superkingdom:Viruses,clade:Riboviria,kingdom:Orthornavirae,phylum:Negarnaviricota,subphylum
287147  superkingdom:Viruses,clade:Riboviria,kingdom:Orthornavirae,phylum:Negarnaviricota,subphylum
***/

static const uint MAX_RECS = 150000;
static const uint MAX_TAXON = 4*1000*1000;

static FILE *g_fOut1;
static FILE *g_fOut2;

static vector<uint> g_Us;
static vector<uint> g_Ss;
static vector<uint> g_Gs;
static vector<uint> g_Fs;
static vector<uint> g_Fams;
static vector<uint> g_Gens;
static vector<uint> g_Spes;

static vector<string> g_SpeNames;
static vector<string> g_GenNames;
static vector<string> g_FamNames;

static map<string, uint> g_SpeToIx;
static map<string, uint> g_GenToIx;
static map<string, uint> g_FamToIx;

static vector<set<uint> > g_SToRecs;
static vector<set<uint> > g_GToRecs;
static vector<set<uint> > g_FToRecs;

static vector<set<uint> > g_SpeToRecs;
static vector<set<uint> > g_GenToRecs;
static vector<set<uint> > g_FamToRecs;

static vector<set<uint> > g_SToSpes;
static vector<set<uint> > g_GToGens;
static vector<set<uint> > g_FToFams;

static vector<set<uint> > g_SpeToSs;
static vector<set<uint> > g_GenToGs;
static vector<set<uint> > g_FamToFs;

static uint g_RecCount;

template<class t> void IncCountMap(map<t, unsigned> &Map, const t &Key, unsigned n = 1)
	{
	if (Map.find(Key) == Map.end())
		Map[Key] = n;
	else
		Map[Key] += n;
	}

template<class t> void CountMapToVecs(const map<t, unsigned> &Map,
  vector<t> &Keys, vector<unsigned> &Counts)
	{
	Keys.clear();
	Counts.clear();
	vector<t> Keys1;
	vector<unsigned> Counts1;
	for (typename map<t, unsigned>::const_iterator p = Map.begin(); p != Map.end(); ++p)
		{
		Keys1.push_back(p->first);
		Counts1.push_back(p->second);
		}
	const unsigned N = SIZE(Keys1);
	if (N == 0)
		return;
	unsigned *Order = myalloc(unsigned, N);
	QuickSortOrderDesc(Counts1.data(), N, Order);
	for (unsigned k = 0; k < N; ++k)
		{
		unsigned i = Order[k];
		Keys.push_back(Keys1[i]);
		Counts.push_back(Counts1[i]);
		}
	myfree(Order);
	}

static uint GetAdd(map<string, uint> &Map, vector<string> &Keys, const string &Key)
	{
	map<string, uint>::iterator p = Map.find(Key);
	uint Ix;
	if (p == Map.end())
		{
		Ix = SIZE(Keys);
		Map[Key] = Ix;
		Keys.push_back(Key);
		}
	else
		Ix = p->second;
	return Ix;
	}

static uint GetIntersectSize(const set<uint> &Set1, const set<uint> &Set2)
	{
	uint Size1 = SIZE(Set1);
	uint Size2 = SIZE(Set2);
	const set<uint> &Big = (Size1 > Size2 ? Set1 : Set2);
	const set<uint> &Small = (Size1 > Size2 ? Set2 : Set1);

	uint Size = 0;
	for (set<uint>::const_iterator pBig = Big.begin(); pBig != Big.end(); ++pBig)
		{
		uint U = *pBig;
		if (Small.find(U) != Small.end())
			++Size;
		}
	return Size;
	}

static void DoS(const string &RankName, uint S,
  const vector<uint> &Spes,
  const set<uint> &SRecs,
  const set<uint> &SSpes,
  const vector<string> &SpeNames,
  const vector<set<uint> > &SpeToRecs)
	{
	const uint RecCount = SIZE(g_Us);
	asserta(SIZE(Spes) == RecCount);

	uint Total = SIZE(SRecs);
	map<uint, uint> SpeToCount;
	uint Total2 = 0;
	uint NamedCladeCount = 0;
	for (set<uint>::const_iterator p = SSpes.begin();
	  p != SSpes.end(); ++p)
		{
		uint Spe = *p;
		asserta(Spe < SIZE(SpeToRecs));
		const string &SpeName = g_SpeNames[Spe];
		if (SpeName != ".")
			++NamedCladeCount;
		const set<uint> &SpeRecs = SpeToRecs[Spe];
		uint n = GetIntersectSize(SRecs, SpeRecs);
		asserta(n > 0);
		Total2 += n;
		asserta(SpeToCount.find(Spe) == SpeToCount.end());
		SpeToCount[Spe] = n;
		}
	bool AllNamed = (NamedCladeCount == SIZE(SSpes));

	char RankChar = toupper(RankName[0]);
	if (Total2 != Total)
		Die("%c%u total2 %u, %u\n", RankChar, S, Total, Total2);

	vector<uint> SpeIxs;
	vector<uint> SpeCounts;
	CountMapToVecs(SpeToCount, SpeIxs, SpeCounts);
	const uint N = SIZE(SpeIxs);
	asserta(SIZE(SpeCounts) == N);

	FILE *f = g_fOut1;
	fprintf(f, "otu=%c%u", RankChar, S);
	fprintf(f, "\tsize=%u", Total);
	fprintf(f, "\tnamed_clades=%u", NamedCladeCount);
	uint Total3 = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint SpeIx = SpeIxs[i];
		uint Count = SpeCounts[i];
		Total3 += Count;
		const char *SpeName = SpeNames[SpeIx].c_str();
		if (strcmp(SpeName, ".") == 0)
			SpeName = "unclassified";
		fprintf(f, "\t%s(%u)", SpeName, Count);
		}
	if (Total3 != Total)
		Die("%c%u total3 %u, %u\n", RankChar, S, Total, Total3);
	fprintf(f, "\n");

	f = g_fOut2;

	string OTUName;
	Ps(OTUName, "%c%u", RankChar, S);
	fprintf(f, "otu=%s", OTUName.c_str());
	fprintf(f, "\totu_size=%u", Total);
	if (SIZE(SpeIxs) == 0)
		fprintf(f, "\tERROR_clades=0");
	else
		{
		asserta(SIZE(SpeCounts) > 0);
		uint TopSpe = SpeIxs[0];
		uint TopSpeCount = SpeCounts[0];
		const string &TopSpeName = SpeNames[TopSpe];
		uint TopSpeSize = SIZE(SpeToRecs[TopSpe]);
		double fOtu = double(TopSpeCount)/Total;
		double fSpe = double(TopSpeCount)/TopSpeSize;
		bool Recip = (fOtu >= 0.5 && fSpe >= 0.5);

		fprintf(f, "\ttop_clade=%s:%s", RankName.c_str(), TopSpeName.c_str());
		fprintf(f, "\ttop_clade_size=%u", TopSpeSize);
		fprintf(f, "\tintersection=%u", TopSpeCount);
		fprintf(f, "\tfract_otu=%.3f", fOtu);
		fprintf(f, "\tfract_clade=%.3f", fSpe);
		fprintf(f, "\tall_named=%c", yon(AllNamed));
		fprintf(f, "\treciprocal=%c", yon(Recip));
		}
	fprintf(f, "\n");
	}

static void DoSpe(const string &RankName,
  uint Spe,
  const string &SpeName,
  const set<uint> &Ss,
  const set<uint> &SpeRecs,
  const vector<set<uint> > &SToRecs)
	{
	const uint SpeSize = SIZE(SpeRecs);

	map<uint, uint> SToCount;
	uint Total2 = 0;
	for (set<uint>::const_iterator p = Ss.begin(); p != Ss.end(); ++p)
		{
		uint S = *p;
		const set<uint> SRecs = SToRecs[S];
		uint n = GetIntersectSize(SpeRecs, SRecs);
		asserta(n > 0);
		Total2 += n;
		asserta(SToCount.find(S) == SToCount.end());
		SToCount[S] = n;
		}
	asserta(Total2 == SpeSize);

	vector<uint> SIxs;
	vector<uint> SCounts;
	CountMapToVecs(SToCount, SIxs, SCounts);
	const uint N = SIZE(SIxs);
	asserta(SIZE(SCounts) == N);
	char RankChar = toupper(RankName[0]);

	FILE *f = g_fOut1;
	fprintf(f, "clade=%s:%s", RankName.c_str(), SpeName.c_str());
	fprintf(f, "\tsize=%u", SpeSize);
	fprintf(f, "\t%c_otus=%u", RankChar, N);
	uint Total3 = 0;
	for (uint i = 0; i < N; ++i)
		{
		uint S = SIxs[i];
		uint Count = SCounts[i];
		Total3 += Count;
		fprintf(f, "\t%c%u(%u)", RankChar, S, Count);
		}
	fprintf(f, "\n");
	asserta(Total3 == SpeSize);

	f = g_fOut2;

	fprintf(f, "clade=%s:%s", RankName.c_str(), SpeName.c_str());
	fprintf(f, "\tclade_size=%u", SpeSize);
	if (SIZE(SIxs) == 0)
		fprintf(f, "\tERROR_otus=0");
	else
		{
		asserta(SIZE(SCounts) > 0);
		uint TopS = SIxs[0];
		uint TopSCount = SCounts[0];
		uint TopSSize = SIZE(SToRecs[TopS]);
		double fOtu = double(TopSCount)/SpeSize;
		double fSpe = double(TopSCount)/TopSSize;
		bool Recip = (fOtu >= 0.5 && fSpe >= 0.5);

		string OTUName;
		Ps(OTUName, "%c%u", RankChar, TopS);

		fprintf(f, "\totu=%s", OTUName.c_str());
		fprintf(f, "\totu_size=%u", TopSSize);
		fprintf(f, "\tintersection=%u", TopSCount);
		fprintf(f, "\tfract_clade=%.3f", fSpe);
		fprintf(f, "\tfract_otu=%.3f", fOtu);
		fprintf(f, "\treciprocal=%c", yon(Recip));
		}
	fprintf(f, "\n");
	}
 
static void Do(const string &RankName,
  const vector<uint> &Ss,
  const vector<uint> &Spes,
  const vector<set<uint> > &SToRecs,
  const vector<set<uint> > &SToSpes,
  const vector<set<uint> > &SpeToSs,
  const vector<string> &SpeNames,
  const vector<set<uint> > &SpeToRecs)
	{
	const uint RecCount = SIZE(Ss);
	const uint SCount = SIZE(SToRecs);
	const uint SpeCount = SIZE(SpeNames);

	set<uint> CheckSet;
	uint RecCount2 = 0;
	ProgressStep(0, SCount, "%s pass 1", RankName.c_str());
	for (uint S = 1; S < SCount; ++S)
		{
		ProgressStep(S, SCount, "%s pass 1", RankName.c_str());
		const set<uint> &Recs = SToRecs[S];
		const uint n = SIZE(Recs);
		if (n == 0)
			continue;
		RecCount2 += n;

		for (set<uint>::const_iterator p = Recs.begin(); p != Recs.end(); ++p)
			{
			uint RecIx = *p;
			asserta(CheckSet.find(RecIx) == CheckSet.end());
			CheckSet.insert(RecIx);
			}

		DoS(RankName, S, Spes, SToRecs[S], SToSpes[S], SpeNames, SpeToRecs);
		}

	for (uint i = 0; i < g_RecCount; ++i)
		if (CheckSet.find(i) == CheckSet.end())
			Warning("%u not found", i);
	if (RecCount2 != RecCount)
		Warning("RecCount2 %s %u, %u (check=%u, g=%u)",
		  RankName.c_str(), RecCount2, RecCount, SIZE(CheckSet), g_RecCount);

	uint RecCount3 = 0;
	for (uint Spe = 0; Spe < SpeCount; ++Spe)
		{
		const set<uint> &Recs = SpeToRecs[Spe];
		RecCount3 += SIZE(Recs);
		ProgressStep(Spe, SpeCount, "%s pass 2", RankName.c_str());

		DoSpe(RankName, Spe, SpeNames[Spe], SpeToSs[Spe], SpeToRecs[Spe], SToRecs);
		}
	if (RecCount3 != RecCount)
		Warning("RecCount3 %u, %u", RecCount2, RecCount);
	}

void OTUMaps()
	{
	const string OTUTsvFN = opt_otumaps;
	const string Out1FN = opt_output;
	const string Out2FN = opt_output2;

	g_fOut1 = CreateStdioFile(Out1FN);
	g_fOut2 = CreateStdioFile(Out2FN);
	asserta(g_fOut1 != 0 && g_fOut2 != 0);

	g_FToRecs.resize(MAX_RECS);
	g_GToRecs.resize(MAX_RECS);
	g_SToRecs.resize(MAX_RECS);

	g_FamToRecs.resize(MAX_RECS);
	g_GenToRecs.resize(MAX_RECS);
	g_SpeToRecs.resize(MAX_RECS);

	g_FamToFs.resize(MAX_RECS);
	g_GenToGs.resize(MAX_RECS);
	g_SpeToSs.resize(MAX_RECS);

	g_FToFams.resize(MAX_RECS);
	g_GToGens.resize(MAX_RECS);
	g_SToSpes.resize(MAX_RECS);

	FILE *f = OpenStdioFile(OTUTsvFN);
	Progress("Reading %s ...", OTUTsvFN.c_str());
	string Line;
	vector<string> Fields;
	ReadLineStdioFile(f, Line);
	asserta(StartsWith(Line, "u	s	g	f"));
	while (ReadLineStdioFile(f, Line))
		{
//      0       1      2       3        4              5           6                   7                 8
// U27330  S10165  G7426   F1290   490041  Caliciviridae   Norovirus       Norwalk virus   Negarnaviricota
		Split(Line, Fields, '\t');
		if (SIZE(Fields) != 10)
			Die("Expected 10 fields in %s", Line.c_str());

		asserta(Fields[0][0] == 'U');
		asserta(Fields[1][0] == 'S');
		asserta(Fields[2][0] == 'G');
		asserta(Fields[3][0] == 'F');

		uint U = (uint) atoi(Fields[0].c_str() + 1);
		uint S = (uint) atoi(Fields[1].c_str() + 1);
		uint G = (uint) atoi(Fields[2].c_str() + 1);
		uint F = (uint) atoi(Fields[3].c_str() + 1);
		uint Taxon = (uint) atoi(Fields[5].c_str());
		if (Taxon == 0 || Taxon > MAX_TAXON)
			Die("Taxon %d in line %s", Taxon, Line.c_str());

		const string &SpeName = Fields[6];
		const string &GenName = Fields[7];
		const string &FamName = Fields[8];

		uint SpeIx = GetAdd(g_SpeToIx, g_SpeNames, SpeName);
		uint GenIx = GetAdd(g_GenToIx, g_GenNames, GenName);
		uint FamIx = GetAdd(g_FamToIx, g_FamNames, FamName);

		uint RecIx = SIZE(g_Us);
		asserta(SIZE(g_Ss) == RecIx);
		asserta(SIZE(g_Gs) == RecIx);
		asserta(SIZE(g_Fs) == RecIx);
		asserta(SIZE(g_Spes) == RecIx);
		asserta(SIZE(g_Gens) == RecIx);
		asserta(SIZE(g_Fams) == RecIx);

		g_Us.push_back(U);
		g_Ss.push_back(S);
		g_Gs.push_back(G);
		g_Fs.push_back(F);

		g_Spes.push_back(SpeIx);
		g_Gens.push_back(GenIx);
		g_Fams.push_back(FamIx);

		asserta(SpeIx < SIZE(g_SpeToRecs));
		asserta(GenIx < SIZE(g_GenToRecs));
		asserta(FamIx < SIZE(g_FamToRecs));

		g_SpeToRecs[SpeIx].insert(RecIx);
		g_GenToRecs[GenIx].insert(RecIx);
		g_FamToRecs[FamIx].insert(RecIx);

		asserta(S < SIZE(g_SToRecs));
		asserta(G < SIZE(g_GToRecs));
		asserta(F < SIZE(g_FToRecs));

		g_SToRecs[S].insert(RecIx);
		g_GToRecs[G].insert(RecIx);
		g_FToRecs[F].insert(RecIx);

		g_SToSpes[S].insert(SpeIx);
		g_GToGens[G].insert(GenIx);
		g_FToFams[F].insert(FamIx);

		g_SpeToSs[SpeIx].insert(S);
		g_GenToGs[GenIx].insert(G);
		g_FamToFs[FamIx].insert(F);
		}
	CloseStdioFile(f);
	Progress("done.\n");

	g_RecCount = SIZE(g_Us);
	asserta(SIZE(g_Ss) == g_RecCount);
	asserta(SIZE(g_Gs) == g_RecCount);
	asserta(SIZE(g_Fs) == g_RecCount);
	asserta(SIZE(g_Spes) == g_RecCount);
	asserta(SIZE(g_Gens) == g_RecCount);
	asserta(SIZE(g_Fams) == g_RecCount);

	Do("family", g_Fs, g_Fams, g_FToRecs, g_FToFams, g_FamToFs, g_FamNames, g_FamToRecs);
	Do("genus", g_Gs, g_Gens, g_GToRecs, g_GToGens, g_GenToGs, g_GenNames, g_GenToRecs);
	Do("species", g_Ss, g_Spes, g_SToRecs, g_SToSpes, g_SpeToSs, g_SpeNames, g_SpeToRecs);

	CloseStdioFile(g_fOut1);
	CloseStdioFile(g_fOut2);
	}
