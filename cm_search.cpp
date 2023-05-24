#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "cmsearcher.h"
#include <time.h>

#if 0

static uint g_DoneCount;
static uint g_HitCount;

static void Thread(ChainReader &CR)
	{
	CMSearcher CS;
	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

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

		const char *Seq = Q.m_Seq.c_str();
		const string &Label = Q.m_Label;
		CS.Search(Q);

		uint APos = UINT_MAX;
		uint BPos = UINT_MAX;
		uint CPos = UINT_MAX;
		double PalmScore = DBL_MAX;
		uint RefIndex = CS.GetPSSMStarts(APos, BPos, CPos, PalmScore);
		if (RefIndex == UINT_MAX)
			continue;
		++g_HitCount;

		if (g_ftsv != 0)
			{
			asserta(RefIndex < SIZE(CS.m_Profs));
			const char *RefLabel = CS.m_Profs[RefIndex]->m_Label.c_str();
			const char *SeqA = Seq + APos;
			const char *SeqB = Seq + BPos;
			const char *SeqC = Seq + CPos;

			char Gate = SeqA[8];
			string GDD;
			GDD += SeqC[2];
			GDD += SeqC[3];
			GDD += SeqC[4];

			fprintf(g_ftsv, "%s", Label.c_str());
			fprintf(g_ftsv, "\t%s", RefLabel);
			fprintf(g_ftsv, "\t%.4f", PalmScore);
			fprintf(g_ftsv, "\t%u", APos+1);
			fprintf(g_ftsv, "\t%u", BPos+1);
			fprintf(g_ftsv, "\t%u", CPos+1);
			fprintf(g_ftsv, "\t%*.*s", AL, AL, SeqA);
			fprintf(g_ftsv, "\t%*.*s", BL, BL, SeqB);
			fprintf(g_ftsv, "\t%*.*s", CL, CL, SeqC);
			fprintf(g_ftsv, "\t%s", (APos < CPos ? "ABC" : "CAB"));
			fprintf(g_ftsv, "\n");
			}
		}
	}

void cmd_cm_search()
	{
	const string &QueryFN = opt_cm_search;

	time_t tStart = time(0);
	if (!optset_db)
		Die("Must specify -db");
	const string &DBFileName = opt_db;

	CMSearcher::ProfsFromFile(DBFileName);

	ChainReader CR;
	CR.Open(QueryFN, false);

	uint ThreadCount = GetRequestedThreadCount();

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR);
	time_t tEnd = time(0);
	uint Secs = uint(tEnd - tStart);
	if (Secs <= 0)
		Secs = 1;
	double Throughput = double(g_DoneCount)/(Secs*ThreadCount);
	ProgressLog("%u done, %u hits, %s secs (%u threads, %.1f/ sec/ thread)\n",
	  g_DoneCount, g_HitCount, IntToStr(Secs), ThreadCount, Throughput);
	}
#else
void cmd_cm_search() { Die("TODO"); }
#endif