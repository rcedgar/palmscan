#include "myutils.h"
#include "outputfiles.h"
#include <map>
#include <set>

static void ParseCat(const string &Cat,
  string &Cls, string &Fold, string &SF, string &Fam)
	{
	vector<string> Fields;
	Split(Cat, Fields, '.');
	asserta(SIZE(Fields) == 4);
	const string &C = Fields[0];
	const string &F = Fields[1];
	const string &S = Fields[2];
	const string &D = Fields[3];

	Cls = C;
	Fold = C + "." + F;
	SF = Fold + "." + S;
	Fam = SF + "." + D;
	}

static int GetFoldTrue(uint Fold1, uint Fold2,
  uint SF1, uint SF2, uint Fam1, uint Fam2)
	{
	if (Fold1 != Fold2)
		return 0;
	if (Fam1 == Fam2)
		{
		asserta(SF1 == SF2);
		return -1;
		}
	if (SF1 == SF2)
		return -1;
	assert(Fold1 == Fold2);
	return 1;
	}

static int GetSFTrue(uint Fold1, uint Fold2,
  uint SF1, uint SF2, uint Fam1, uint Fam2)
	{
	if (Fold1 != Fold2)
		return 0;
	if (Fam1 == Fam2)
		return -1;
	else if (SF1 == SF2)
		return 1;
	else if (Fold1 == Fold2 && SF1 != SF2)
		return -1;
	asserta(false);
	return -9;
	}

static int GetFamTrue(uint Fold1, uint Fold2,
  uint SF1, uint SF2, uint Fam1, uint Fam2)
	{
	if (Fold1 != Fold2)
		return 0;
	else if (Fam1 == Fam2)
		return 1;
	else
		return -1;
	}

static uint g_ScoreType;
static uint ScoreToBin(double Score)
	{
	switch (g_ScoreType)
		{
	case 0: // TMalign, 0.0 to 1.0
		{
		asserta(Score >= 0 && Score <= 1);
		double r = Score*100.0;
		uint Bin = uint(r);
		asserta(Bin >= 0 && Bin <= 100);
		return Bin;
		}

	case 1: // foldseek
		{
		if (Score < 1e-40)
			Score = 1e-40;
		if (Score > 10)
			Score = 10;

	// range -300 to +2
		double LogScore = log10(Score);
		double r = -(LogScore - 2.0)*100.0/42.0;
		if (r < -0.1 ||  r > 100.1)
			Die("Score = %.4g, LogScore = %.4g, r = %.4g",
			  Score, LogScore, r);
		if (r < 0)
			r = 0;
		uint Bin = uint(r);
		if (Bin > 100)
			Bin = 100;
		return Bin;
		}

	case 2: // DALI
		{
		asserta(Score >= 0 && Score <= 90);
		if (Score > 20)
			Score = 20;
		uint Bin = uint(Score*5.0);
		asserta(Bin >= 0 && Bin <= 100);
		return Bin;
		}

	case 3: // 3dblast
		{
		asserta(Score >= 50 && Score <= 10000);
		Score -= 50;
		asserta(Score >= 0);
		if (Score > 200)
			Score = 200;
		uint Bin = uint(Score/2);
		asserta(Bin >= 0 && Bin <= 100);
		return Bin;
		}

	case 4: // 3dblastsw
		{
		if (Score < 0)
			Score = 0;
		asserta(Score >= 0);
		if (Score > 250)
			Score = 250;
		uint Bin = uint(Score/2.5);
		asserta(Bin >= 0 && Bin <= 100);
		return Bin;
		}

	case 5: // ce
		{
		asserta(Score >= 0 && Score <= 10);
		uint Bin = uint(Score*10);
		asserta(Bin >= 0 && Bin <= 100);
		return Bin;
		}

		} // end switch

	asserta(false);
	return UINT_MAX;
	}

void cmd_scop40_firstfp()
	{
	const string &IdClassFN = opt_scop40_firstfp;
	const string &AlnTsvFN = opt_input;

	g_ScoreType = opt_scoretype;
	uint ScoreFieldNr = 2;
	if (opt_scoretype == 1)
		ScoreFieldNr = 10;

	vector<string> Domains;
	vector<string> Cats;
	map<string, uint > CatToIndex;
	map<string, uint > DomainToIndex;
	vector<uint> DomainIndexToCatIndex;
	map<uint, vector<uint> > CatToDomains;

	vector<string> Folds;
	map<string, uint> FoldToIndex;

	vector<string> SFs;
	map<string, uint> SFToIndex;

	vector<string> Fams;
	map<string, uint> FamToIndex;

	string Line;
	vector<string> Fields;

	vector<uint> DomainToFold;
	vector<uint> DomainToSF;
	vector<uint> DomainToFam;

	FILE *f = OpenStdioFile(IdClassFN);
	Progress("Reading... ");
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2);

		const string &Domain = Fields[0];
		const string &Cat = Fields[1];

		string Cls, Fold, SF, Fam;
		ParseCat(Cat, Cls, Fold, SF, Fam);

		uint FoldIndex = UINT_MAX;
		if (FoldToIndex.find(Fold) == FoldToIndex.end())
			{
			FoldIndex = SIZE(Folds);
			Folds.push_back(Fold);
			FoldToIndex[Fold] = FoldIndex;
			}
		else
			FoldIndex = FoldToIndex[Fold];

		uint SFIndex = UINT_MAX;
		if (SFToIndex.find(SF) == SFToIndex.end())
			{
			SFIndex = SIZE(SFs);
			SFs.push_back(SF);
			SFToIndex[SF] = SFIndex;
			}
		else
			SFIndex = SFToIndex[SF];

		uint FamIndex = UINT_MAX;
		if (FamToIndex.find(Fam) == FamToIndex.end())
			{
			FamIndex = SIZE(Fams);
			Fams.push_back(Fam);
			FamToIndex[Fam] = FamIndex;
			}
		else
			FamIndex = FamToIndex[Fam];

		uint CatIndex = UINT_MAX;
		map<string, uint>::const_iterator p = CatToIndex.find(Cat);
		if (p == CatToIndex.end())
			{
			CatIndex = SIZE(Cats);
			CatToIndex[Cat] = CatIndex;
			Cats.push_back(Cat);
			CatToDomains[CatIndex].resize(CatIndex+1);
			}
		else
			CatIndex = p->second;

		if (DomainToIndex.find(Domain) != DomainToIndex.end())
			Die("Duplicated domain '%s'", Domain.c_str());

		uint DomainIndex = SIZE(Domains);

		Domains.push_back(Domain);
		DomainToIndex[Domain] = DomainIndex;
		DomainIndexToCatIndex.push_back(CatIndex);
		CatToDomains[CatIndex].push_back(DomainIndex);

		DomainToFold.push_back(FoldIndex);
		DomainToSF.push_back(SFIndex);
		DomainToFam.push_back(FamIndex);
		}
	CloseStdioFile(f);
	Progress("done.\n");

	const uint DomainCount = SIZE(Domains);
	asserta(SIZE(DomainToFold) == DomainCount);
	asserta(SIZE(DomainToSF) == DomainCount);
	asserta(SIZE(DomainToFam) == DomainCount);
	asserta(SIZE(DomainIndexToCatIndex) == DomainCount);

	f = OpenStdioFile(AlnTsvFN);
	uint64 FileSize = GetStdioFileSize64(f);
	uint LineNr = 0;

	char Delim = ' ';
	if (g_ScoreType == 1 || g_ScoreType == 6)
		Delim = '\t';
	uint CurrentDom = UINT_MAX;
	uint CurrentDomFPCount = 0;
	uint CurrentDomTPCount = 0;
	uint SumAllTPCount = 0;
	uint SumFirstTPCount = 0;
	double LastScore = DBL_MAX;
	while (ReadLineStdioFile(f, Line))
		{
		++LineNr;
		if (LineNr%100000 == 0)
			{
			uint64 FilePos = GetStdioFilePos64(f);
			double Pct = FilePos*100.0/FileSize;
			Progress("%.1f%%  %s\r", Pct, AlnTsvFN.c_str());
			}
		if (SIZE(Line) == 0)
			continue;
		if (Line[0] == ' ')
			continue;
		Split(Line, Fields, Delim);
		uint n = SIZE(Fields);
		if (ScoreFieldNr >= n)
			continue;
		asserta(ScoreFieldNr < n);
		const string &Domain1 = Fields[0];
		const string &Domain2 = Fields[1];
		if (Domain1 == Domain2)
			continue;

		double Score = StrToFloat(Fields[ScoreFieldNr]);

		map<string, uint>::const_iterator p1 = DomainToIndex.find(Domain1);
		map<string, uint>::const_iterator p2 = DomainToIndex.find(Domain2);

		if (p1 == DomainToIndex.end())
			continue;
		if (p2 == DomainToIndex.end())
			continue;

		uint Dom1 = p1->second;
		uint Dom2 = p2->second;

		uint Fold1 = DomainToFold[Dom1];
		uint Fold2 = DomainToFold[Dom2];
		bool TP = (Fold1 == Fold2);

		if (g_ftsv != 0)
			{
			FILE *f = g_ftsv;

			uint Fam1 = DomainToFam[Dom1];
			uint Fam2 = DomainToFam[Dom2];

			fprintf(f, "%s", Domain1.c_str());
			fprintf(f, "\t%s", Domain2.c_str());
			fprintf(f, "\t%s", Fams[Fam1].c_str());
			fprintf(f, "\t%s", Fams[Fam2].c_str());
			fprintf(f, "\t%.4g", Score);
			fprintf(f, "\t%u", CurrentDomFPCount);
			fprintf(f, "\t%c", tof(TP));
			fprintf(f, "\n");
			}

		if (Dom1 != CurrentDom)
			{
			CurrentDom = Dom1;
			CurrentDomTPCount = 0;
			CurrentDomFPCount = 0;
			}
		else
			asserta(Score <= LastScore);

		if (TP)
			++SumAllTPCount;
		else
			{
			if (CurrentDomFPCount == 0)
				++SumFirstTPCount;
			++CurrentDomFPCount;
			}
		}
	Log("Sum all %u, first %u\n", SumAllTPCount, SumFirstTPCount);
	}
