#include "myutils.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "cddata.h"
#include "cdsearcher.h"
#include "chainreader.h"
#include "cmp.h"
#include "omplock.h"

static uint g_DoneCount;
static uint g_HitCount;

void SetPalmTemplate(const CDInfo &Info, CDTemplate &Tpl);

static void HitToTsv(const CDSearcher &CS, double Score,
  const vector<uint> &MotifCoords)
	{
	FILE *f = g_ftsv;
	if (f == 0)
		return;
	if (Score == 0)
		return;
	if (MotifCoords.empty())
		return;

	const PDBChain &Q = *CS.m_Query;
	const char *Seq = Q.m_Seq.c_str();

	double P_RdRp = 0;
	uint PosA = MotifCoords[1];
	uint PosC = MotifCoords[3];
	char Gate = '.';
	string GDD = "...";
	if (PosA != UINT_MAX)
		{
		asserta(PosA + 8 < SIZE(Q.m_Seq));
		Gate = Seq[PosA + 8];
		}

	if (PosC != UINT_MAX)
		{
		asserta(PosC + 4 < SIZE(Q.m_Seq));
		GDD[0] = Seq[PosC + 2];
		GDD[1] = Seq[PosC + 3];
		GDD[2] = Seq[PosC + 4];
		}

	if (PosA != UINT_MAX && PosC != UINT_MAX)
		P_RdRp = CMP::GetRdRpProb(Gate, GDD);

	const uint MotifCount = CS.m_Info->GetMotifCount();
	asserta(SIZE(MotifCoords) == MotifCount);
	vector<uint> MotifIndexes;

	double FinalScore = P_RdRp*Score;
	if (FinalScore < 0.1)
		return;

	Lock();
	fprintf(f, "%.4f", FinalScore);
	fprintf(f, "\t%.4f", Score);
	fprintf(f, "\t%.4f", P_RdRp);
	fprintf(f, "\t%s", Q.m_Label.c_str());
	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		{
		MotifIndexes.push_back(MotifIndex);
		uint ML = CS.m_Info->GetMotifLength(MotifIndex);
		uint Pos = MotifCoords[MotifIndex];
		if (Pos == UINT_MAX)
			{
			fprintf(f, "\t.\t.");
			continue;
			}

		fprintf(f, "\t%u\t%*.*s", Pos + 1, ML, ML, Seq + Pos);
		}
	fprintf(f, "\t%c\t%s", Gate, GDD.c_str());
	fprintf(f, "\n");
	Unlock();
	}

static void Thread(ChainReader &CR, const CDInfo &Info, 
  const CDTemplate &Tpl, const CDData &DataAvg,
  const CDData &DataStdDev)
	{
	CDSearcher CS;
	CS.Init(Info, DataAvg, DataStdDev);
	CS.m_Template = &Tpl;

	PDBChain Q;
	vector<uint> Hit;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		double Score = CS.SearchPalm(Q, Hit);
		if (Score > 0.5)
			++g_HitCount;
		HitToTsv(CS, Score, Hit);

#pragma omp critical
		{
		if (++g_DoneCount%1000 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}
		}
		}
	}

void cmd_cdp_search()
	{
	const string &InputFN = opt_cdp_search;

	asserta(optset_model);
	FILE *f = OpenStdioFile(opt_model);
	CDInfo Info;
	Info.FromTsv(f);

	CDData DataAvg;
	CDData DataStdDev;
	DataAvg.FromTsv(Info, f);
	DataStdDev.FromTsv(Info, f);

	CDTemplate Tpl;
	SetPalmTemplate(Info, Tpl);

	ChainReader CR;
	CR.Open(InputFN);
	uint ThreadCount = GetRequestedThreadCount();

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Info, Tpl, DataAvg, DataStdDev);

	ProgressLog("100.0%% done, %u / %u hits\n",
	  g_HitCount, g_DoneCount);
	}
