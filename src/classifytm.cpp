#include "myutils.h"
#include "pdbchain.h"
#include "calreader.h"
#include "outputfiles.h"
#include "omplock.h"
#include "tma.h"

static const double STRONG_RDRP = 0.600;
static const double WEAK_RDRP = 0.550;
static const double POSSIBLE_RDRP = 0.500;
static const double WEAK_DIFF = 0.1;

static uint g_RecCount;
static uint g_RdRpCount;
static uint g_OtherCount;

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

static void WriteAccAndTM(FILE *f, uint Index,
  vector<PDBChain *> &Rs, double TM)
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

	fprintf(f, "\t%s\t%.3f", Acc.c_str(), TM);
	}

static void Classify1(FILE *f, PDBChain &Q, vector<PDBChain *> &Rs,
  const vector<bool> &RefIsRdRp, TMA &T)
	{
	const uint NR = SIZE(Rs);
	double TM_TopRdRp = DBL_MAX;
	double TM_TopOther = DBL_MAX;
	uint RefIndex_TopRdRp = UINT_MAX;
	uint RefIndex_TopOther = UINT_MAX;
	vector<double> TMs;
	double TopTM = 0;
	uint TopIndex = UINT_MAX;
	for (uint RefIndex = 0; RefIndex < NR; ++RefIndex)
		{
		const PDBChain &R = *Rs[RefIndex];
		double TM = T.AlignChains(Q, R);
		TMs.push_back(TM);
		if (TM > TopTM)
			{
			TopTM = TM;
			TopIndex = RefIndex;
			}
		}
	asserta(TopIndex != UINT_MAX);
	bool TopIsRdRp = RefIsRdRp[TopIndex];
	double SecondTM = 999;
	uint SecondIndex = UINT_MAX;
	for (uint RefIndex = 0; RefIndex < NR; ++RefIndex)
		{
		if (RefIsRdRp[RefIndex] != TopIsRdRp &&
		  TMs[RefIndex] < SecondTM)
			{
			SecondIndex = RefIndex;
			SecondTM = TMs[RefIndex];
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

	double PalmScore = TopTM;
	double RdRpScore = 0;
	if (TopIsRdRp)
		{
		RdRpScore = PalmScore;
		if (SecondTM != DBL_MAX)
			{
			double SecondScore = SecondTM;
			DampScore(RdRpScore, SecondScore/2);
			}
		}

	if (!TopIsRdRp)
		{
		double TMDiff = SecondTM - TopTM;
		if (TMDiff <= WEAK_DIFF)
			{
			RdRpScore = SecondTM;
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

	string Cat = "other";
	if (RdRpScore >= STRONG_RDRP)
		Cat = "rdrp.strong";
	else if (RdRpScore >= WEAK_RDRP)
		Cat = "rdrp.weak";
	else if (RdRpScore >= POSSIBLE_RDRP)
		Cat = "rdrp.possible";

	uint BPos = Q.m_MotifPosVec[1];
	const uint PPL = Q.GetSeqLength();

	Lock();
	static bool HdrDone = false;
	if (!HdrDone)
		{
		fprintf(f, "Srdrp");
		fprintf(f, "\tClass");
		fprintf(f, "\tTMpalm");
		fprintf(f, "\tQuery");
		fprintf(f, "\tA");
		fprintf(f, "\tB");
		fprintf(f, "\tC");
		fprintf(f, "\tDGD");
		fprintf(f, "\tGDD");
		fprintf(f, "\tTopHit");
		fprintf(f, "\tTMtop");
		fprintf(f, "\tHit2");
		fprintf(f, "\tTM2");
		fprintf(f, "\tPPL");
		fprintf(f, "\tPosB");
		fprintf(f, "\tdAB");
		fprintf(f, "\tdBC");
		fprintf(f, "\tdAC");
		fputc('\n', f);

		HdrDone = true;
		}

	fprintf(f, "%.3f", RdRpScore);
	fprintf(f, "\t%s", Cat.c_str());
	fprintf(f, "\t%.3f", PalmScore);
	WriteSplitLabel(f, &Q);
	fprintf(f, "\t%s", SuperMotif.c_str());
	fprintf(f, "\t%s", GDD.c_str());
	WriteAccAndTM(f, TopIndex, Rs, TopTM);
	WriteAccAndTM(f, SecondIndex, Rs, SecondTM);
	fprintf(f, "\t%u", PPL);
	fprintf(f, "\t%u", BPos);
	fprintf(f, "\t%.1f", dAB);
	fprintf(f, "\t%.1f", dBC);
	fprintf(f, "\t%.1f", dAC);
	fputc('\n', f);
	Unlock();
	}

static void Thread(CalReader &CR, vector<PDBChain *> &Qs,
  vector<PDBChain *> &Rs, const vector<bool> &RefIsRdRp)
	{
	const uint ThreadIndex = GetThreadIndex();
	PDBChain Q;
//	PPCAligner PA;
	TMA T;
	for (;;)
		{
		Lock();
		bool Ok = CR.GetNext(Q);
		Unlock();
		if (!Ok)
			break;
		if (g_RecCount%100 == 0)
			{
			Lock();
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u RdRp hits\r",
			  sPct.c_str(), g_RdRpCount, g_RecCount);
			Unlock();
			}
		Lock();
		++g_RecCount;
		Unlock();
		if (Q.m_MotifPosVec.empty())
			{
			Warning("Missing motif coords >%s", Q.m_Label.c_str());
			continue;
			}
		//PA.SetQuery(Q);
		Classify1(g_ftsv, Q, Rs, RefIsRdRp, T);
		}
	}

void cmd_classifytm()
	{
	const string &QFN = opt_classifytm;
	const string &RefFN = opt_ref;

	vector<PDBChain *> Rs;
	ReadChains(RefFN, Rs);
	const uint NR = SIZE(Rs);
	asserta(NR > 0);

	vector<bool> RefIsRdRp;
	Log("Refs:\n");
	for (uint i = 0; i < NR; ++i)
		{
		const string &Label = Rs[i]->m_Label;
		bool IsRdRp = StartsWith(Label, "rdrp.");
		RefIsRdRp.push_back(IsRdRp);
		Log("%c  %s\n", pom(IsRdRp), Label.c_str());
		}

	uint ThreadCount = GetRequestedThreadCount();

	CalReader CR;
	CR.Open(QFN);
	vector<PDBChain *> Qs(ThreadCount);

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Qs, Rs, RefIsRdRp);

	Progress("100%% done, %u / %u RdRp hits\n",
		g_RdRpCount, g_RecCount);

	Log("%u / %u RdRp hits (%.1f%%)\n",
	  g_RdRpCount, g_RecCount, GetPct(g_RdRpCount, g_RecCount));
	}
