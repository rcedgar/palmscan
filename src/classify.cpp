#include "myutils.h"
#include "pdbchain.h"
#include "ppcaligner.h"
#include "calreader.h"
#include "outputfiles.h"

static const double WEAK_DIFF = 0.5;
static const double RMSD_Score_THRESHOLD = 4.0;

static uint g_RdRpCount;
static uint g_OtherCount;

/***
three_point_quadratic.py
  (0, 1)
  (4, 0.5)
  (6, 0.25)

         x    Score
  --------  ------
         0   1.0000
       0.5   0.9375
         1   0.8750
       1.5   0.8125
         2   0.7500
       2.5   0.6875
         3   0.6250
       3.5   0.5625
         4   0.5000
       4.5   0.4375
         5   0.3750
       5.5   0.3125
         6   0.2500
       6.5   0.1875
         7   0.1250
       7.5   0.0625
         8   0.0000
***/
static double GetRMSDScore(double x)
	{
	static const double a = 0;
	static const double b = -0.125;
	static const double c = 1;
	double Score = a*x*x + b*x + c;
	if (Score < 0)
		Score = 0;
	else if (Score > 1)
		Score = 1;
	return Score;
	}

static double BoostScore(double Score, double x)
	{
	asserta(Score >= 0 && Score <= 1);
	asserta(x >= 0 && x <= 1);
	double d = 1.0 - Score;
	double BoostedScore = Score + d*x;
	asserta(BoostedScore >= 0 && BoostedScore < 1.001);
	if (BoostedScore > 1)
		BoostedScore = 1;
	return BoostedScore;
	}

static double DampScore(double Score, double x)
	{
	asserta(Score >= 0 && Score <= 1);
	asserta(x >= 0 && x <= 1);
	double DampedScore = Score*(1 - x);
	asserta(DampedScore >= 0 && DampedScore <= 1);
	return DampedScore;
	}

static void WriteSplitLabel(FILE *f, const PDBChain *Chain)
	{
	if (f == 0)
		return;

	if (Chain == 0)
		fprintf(f, "\t.\t.\t.\t.");
	else
		{
		string Label = Chain->m_Label;
		vector<string> Fields;
		Split(Label, Fields, ' ');
		asserta(SIZE(Fields) == 4);
		string Acc = Fields[0];
		size_t n = Acc.find("_unrelaxed_rank");
		if (n != string::npos)
			Acc = Acc.substr(0, n);
		fprintf(f, "\t%s", Acc.c_str());

		string ASeq;
		string BSeq;
		string CSeq;
		Chain->GetMotifSeq(0, ASeq);
		Chain->GetMotifSeq(1, BSeq);
		Chain->GetMotifSeq(2, CSeq);

		fputc('\t', f);
		fputs(ASeq.c_str(), f);
		fputc('\t', f);
		fputs(BSeq.c_str(), f);
		fputc('\t', f);
		fputs(CSeq.c_str(), f);
		}
	}

static void WriteAccAndRMSD(FILE *f, uint Index,
  vector<PDBChain *> &Rs, double RMSD)
	{
	if (f == 0)
		return;

	if (Index == UINT_MAX)
		{
		fprintf(f, "\t.\t.");
		return;
		}

	asserta(Index < SIZE(Rs));
	const string &Label = Rs[Index]->m_Label;
	size_t n = Label.find(' ');
	string Acc;
	if (n == string::npos)
		Acc = Label;
	else
		Acc = Label.substr(0, n);

	fprintf(f, "\t%s\t%.3g", Acc.c_str(), RMSD);
	}

static void WriteCat(FILE *f, bool TopIsRdRp,
  double TopRMSD, double SecondRMSD)
	{
	if (f == 0)
		return;

	const char *Cat = ".";
	if (TopRMSD == DBL_MAX)
		Cat = "other.nohit";
	else if (SecondRMSD == DBL_MAX)
		{
		if (TopIsRdRp)
			Cat = "rdrp.no-other-hit";
		else
			Cat = "other.no-rdrp-hit";
		}
	else
		{
		asserta(TopRMSD != DBL_MAX && SecondRMSD != DBL_MAX);
		bool Weak = (SecondRMSD - TopRMSD <= WEAK_DIFF);
		if (TopIsRdRp)
			Cat = "rdrp.top-hit";
		else
			{
			if (Weak)
				Cat = "rdrp-like";
			else
				Cat = "other.top-hit";
			}
		}
	fprintf(f, "\t%s", Cat);
	}

static void Classify1(FILE *f, PDBChain &Q, vector<PDBChain *> &Rs,
  const vector<bool> &RefIsRdRp, PPCAligner &PA)
	{
	const uint NR = SIZE(Rs);
	double RMSD_TopRdRp = DBL_MAX;
	double RMSD_TopOther = DBL_MAX;
	uint RefIndex_TopRdRp = UINT_MAX;
	uint RefIndex_TopOther = UINT_MAX;
	vector<double> RMSDs;
	double TopRMSD = 999;
	uint TopIndex = UINT_MAX;
	for (uint RefIndex = 0; RefIndex < NR; ++RefIndex)
		{
		const PDBChain &R = *Rs[RefIndex];
		PA.SetRef(R);
		double RMSD = PA.GetMotifRMSD();
		RMSDs.push_back(RMSD);
		if (RMSD < TopRMSD)
			{
			TopRMSD = RMSD;
			TopIndex = RefIndex;
			}
		}
	asserta(TopIndex != UINT_MAX);
	bool TopIsRdRp = RefIsRdRp[TopIndex];
	double SecondRMSD = 999;
	uint SecondIndex = UINT_MAX;
	for (uint RefIndex = 0; RefIndex < NR; ++RefIndex)
		{
		if (RefIsRdRp[RefIndex] != TopIsRdRp &&
		  RMSDs[RefIndex] < SecondRMSD)
			{
			SecondIndex = RefIndex;
			SecondRMSD = RMSDs[RefIndex];
			}
		}

	if (TopIsRdRp)
		++g_RdRpCount;
	else
		++g_OtherCount;

	if (f == 0)
		return;
	if (TopIndex == UINT_MAX)
		return;

	double dAB, dBC, dAC;
	Q.GetMotifDists(dAB, dBC, dAC);

	string SuperMotif;
	Q.GetSuperMotif(SuperMotif);

	string GDD;
	Q.GetBestMatchGDD(GDD);

	double PalmScore = GetRMSDScore(TopRMSD);
	double RdRpScore = 0;
	if (TopIsRdRp)
		{
		RdRpScore = PalmScore;
		if (SecondRMSD != DBL_MAX)
			{
			double SecondScore = GetRMSDScore(SecondRMSD);
			DampScore(RdRpScore, SecondScore/2);
			}
		}

	if (!TopIsRdRp)
		{
		double RMSDDiff = SecondRMSD - TopRMSD;
		if (RMSDDiff <= WEAK_DIFF)
			{
			RdRpScore = GetRMSDScore(SecondRMSD);
			RdRpScore = DampScore(RdRpScore, 0.5);
			RdRpScore = DampScore(RdRpScore, PalmScore/2);
			}
		}

	if (GDD == "ADD" || GDD == "LDD" || GDD == "VDD")
		RdRpScore = 0.0;
	else if (GDD == "GDD" || GDD == "SDD" || GDD == "GDN")
		{
		if (SuperMotif == "DGD")
			RdRpScore = BoostScore(RdRpScore, 0.5);
		}

	uint BPos = Q.m_MotifPosVec[1];
	const uint PPL = Q.GetSeqLength();

	fprintf(f, "%.3f", RdRpScore);
	fprintf(f, "\t%.3f", PalmScore);
	WriteSplitLabel(f, &Q);
	fprintf(f, "\t%s", SuperMotif.c_str());
	fprintf(f, "\t%s", GDD.c_str());
	WriteAccAndRMSD(f, TopIndex, Rs, TopRMSD);
	WriteAccAndRMSD(f, SecondIndex, Rs, SecondRMSD);
	fprintf(f, "\t%u", PPL);
	fprintf(f, "\t%u", BPos);
	fprintf(f, "\t%.1f", dAB);
	fprintf(f, "\t%.1f", dBC);
	fprintf(f, "\t%.1f", dAC);
	fputc('\n', f);
	}

void cmd_classify()
	{
	const string &QFN = opt_classify;
	const string &RefFN = opt_ref;

	vector<PDBChain *> Rs;
	ReadChains(RefFN, Rs);
	const uint NR = SIZE(Rs);
	asserta(NR > 0);

	vector<bool> RefIsRdRp;
	for (uint i = 0; i < NR; ++i)
		{
		const string &Label = Rs[i]->m_Label;
		bool IsRdRp = StartsWith(Label, "rdrp.");
		RefIsRdRp.push_back(IsRdRp);
		Log("%c  %s\n", pom(IsRdRp), Label.c_str());
		}

	PPCAligner PA;

	CalReader CR;
	CR.Open(QFN);

	uint RecCount = 0;
	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;
		if (RecCount%10000 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u RdRp hits\r",
			  sPct.c_str(), g_RdRpCount, RecCount);
			}
		++RecCount;
		if (Q.m_MotifPosVec.empty())
			{
			Warning("Missing motif coords >%s", Q.m_Label.c_str());
			continue;
			}
		PA.SetQuery(Q);
		Classify1(g_ftsv, Q, Rs, RefIsRdRp, PA);
		}
	Progress("100%% done, %u / %u RdRp hits\n",
		g_RdRpCount, RecCount);
	}
