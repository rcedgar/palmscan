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

void cmd_scop40()
	{
	const string &IdClassFN = opt_scop40;
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

	if (opt_calctps)
		{
		uint PairCount = 0;
		uint NFold = 0;
		uint NSF = 0;
		uint NFam = 0;

		uint NFoldTrue = 0;
		uint NSFTrue = 0;
		uint NFamTrue = 0;

		for (uint i = 0; i < DomainCount; ++i)
			{
			uint Fold_i = DomainToFold[i];
			uint SF_i = DomainToSF[i];
			uint Fam_i = DomainToFam[i];

			for (uint j = 0; j < DomainCount; ++j)
				{
				if (i == j)
					continue;

				++PairCount;

				uint Fold_j = DomainToFold[j];
				uint SF_j = DomainToSF[j];
				uint Fam_j = DomainToFam[j];

				int FoldTrue = GetFoldTrue(Fold_i, Fold_j, SF_i, SF_j, Fam_i, Fam_j);
				int SFTrue = GetSFTrue(Fold_i, Fold_j, SF_i, SF_j, Fam_i, Fam_j);
				int FamTrue = GetFamTrue(Fold_i, Fold_j, SF_i, SF_j, Fam_i, Fam_j);

				if (FoldTrue == 0 || FoldTrue == 1)
					{
					++NFold;
					if (FoldTrue == 1)
						++NFoldTrue;
					}

				if (SFTrue == 0 || SFTrue == 1)
					{
					++NSF;
					if (SFTrue == 1)
						++NSFTrue;
					}

				if (FamTrue == 0 || FamTrue == 1)
					{
					++NFam;
					if (FamTrue == 1)
						++NFamTrue;
					}
				}
			}

		ProgressLog("Fold %u / %u\n", NFoldTrue, NFold);
		ProgressLog("SF %u / %u\n", NSFTrue, NSF);
		ProgressLog("Fam %u / %u\n", NFamTrue, NFam);
		return;
		}

	if (g_freport_3d != 0)
		{
		FILE *f = g_freport_3d;
		for (uint i = 0; i < SIZE(Folds); ++i)
			fprintf(f, "@fold\t%u\t%s\n",
			  i, Folds[i].c_str());

		for (uint i = 0; i < SIZE(SFs); ++i)
			fprintf(f, "@sf\t%u\t%s\n",
			  i, SFs[i].c_str());

		for (uint i = 0; i < SIZE(Fams); ++i)
			fprintf(f, "@fam\t%u\t%s\n",
			  i, Fams[i].c_str());

		for (uint DomainIndex = 0; DomainIndex < SIZE(Domains); ++DomainIndex)
			{
			fprintf(f, "@domain\t%u\t%s",
			  DomainIndex, Domains[DomainIndex].c_str());

			uint CatIndex = DomainIndexToCatIndex[DomainIndex];
			const string &Cat = Cats[CatIndex];

			fprintf(f, "\t%u\t%s", CatIndex, Cat.c_str());

			uint Fold = DomainToFold[DomainIndex];
			uint SF = DomainToSF[DomainIndex];
			uint Fam = DomainToFam[DomainIndex];

			fprintf(f, "\t%u\t%s", Fold, Folds[Fold].c_str());
			fprintf(f, "\t%u\t%s", SF, SFs[SF].c_str());
			fprintf(f, "\t%u\t%s", Fam, Fams[Fam].c_str());
			fprintf(f, "\n");
			}
		}

	f = OpenStdioFile(AlnTsvFN);
	uint64 FileSize = GetStdioFileSize64(f);
	uint LineNr = 0;
	set<string> MissingDomains;
	uint FoundPairCount = 0;
	uint MissingPairCount = 0;
	vector<uint> Dom1s;
	vector<uint> Dom2s;
	vector<uint> Cat1s;
	vector<uint> Cat2s;

	vector<uint> BinToFoldTrues(101);
	vector<uint> BinToFoldFalses(101);

	vector<uint> BinToSFTrues(101);
	vector<uint> BinToSFFalses(101);

	vector<uint> BinToFamTrues(101);
	vector<uint> BinToFamFalses(101);

	char Delim = ' ';
	if (g_ScoreType == 1)
		Delim = '\t';
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
			MissingDomains.insert(Domain1);
		if (p2 == DomainToIndex.end())
			MissingDomains.insert(Domain2);
		if (p1 == DomainToIndex.end() || p2 == DomainToIndex.end())
			{
			++MissingPairCount;
			continue;
			}
		++FoundPairCount;

		uint Dom1 = p1->second;
		uint Dom2 = p2->second;

		uint Fold1 = DomainToFold[Dom1];
		uint Fold2 = DomainToFold[Dom2];

		uint SF1 = DomainToSF[Dom1];
		uint SF2 = DomainToSF[Dom2];

		uint Fam1 = DomainToFam[Dom1];
		uint Fam2 = DomainToFam[Dom2];

		int FoldTrue = GetFoldTrue(Fold1, Fold2, SF1, SF2, Fam1, Fam2);
		int SFTrue = GetSFTrue(Fold1, Fold2, SF1, SF2, Fam1, Fam2);
		int FamTrue = GetFamTrue(Fold1, Fold2, SF1, SF2, Fam1, Fam2);

		uint Bin = ScoreToBin(Score);

		if (FoldTrue == 1)
			BinToFoldTrues[Bin] += 1;
		else if (FoldTrue == 0)
			BinToFoldFalses[Bin] += 1;

		if (SFTrue == 1)
			BinToSFTrues[Bin] += 1;
		else if (SFTrue == 0)
			BinToSFFalses[Bin] += 1;

		if (FamTrue == 1)
			BinToFamTrues[Bin] += 1;
		else if (FamTrue == 0)
			BinToFamFalses[Bin] += 1;

#if 0
		{
		const string &sDom1 = Domains[Dom1];
		const string &sDom2 = Domains[Dom2];

		const string &sFam1 = Fams[Fam1];
		const string &sFam2 = Fams[Fam2];
		static uint n;
		if (++n > 100)
			break;
		Log("%12.12s  %12.12s  %+3d  %+3d  %+3d\n",
		  sFam1.c_str(), sFam2.c_str(), FoldTrue, SFTrue, FamTrue);
		}
#endif
		}
	Progress("100.0%% %s\n", AlnTsvFN.c_str());
	Progress("%u found pairs, %u missing pairs\n",
	  FoundPairCount, MissingPairCount);

	if (g_ftsv != 0)
		{
		for (uint Bin = 0; Bin <= 100; ++Bin)
			{
			fprintf(g_ftsv, "%.4f", Bin/100.0);
			fprintf(g_ftsv, "\t%u", BinToFoldTrues[Bin]);
			fprintf(g_ftsv, "\t%u", BinToFoldFalses[Bin]);
			fprintf(g_ftsv, "\t%u", BinToSFTrues[Bin]);
			fprintf(g_ftsv, "\t%u", BinToSFFalses[Bin]);
			fprintf(g_ftsv, "\t%u", BinToFamTrues[Bin]);
			fprintf(g_ftsv, "\t%u", BinToFamFalses[Bin]);
			fprintf(g_ftsv, "\n");
			}
		}
	}
