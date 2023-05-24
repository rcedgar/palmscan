#include "myutils.h"
#include "pssm.h"
#include "rdrpsearcher.h"
#include "seqinfo.h"
#include "fastaseqsource.h"
#include "omplock.h"
#include <time.h>

/***
***********************************
TOO MANY FPs
***********************************
***/

static uint g_QueryCount;
static uint g_FoundCount;

void SetExcludes();

static void SearchC(const string &QueryFileName)
	{
	asserta(optset_model);
	bool Trace = opt_trace;
	double MinPalmScore = 10.0;
	double MinCScore = 2.0;
	if (optset_mincscore)
		MinCScore = opt_minscore;
	const string &ModelFileName = opt_model;
	if (optset_min_palm_score)
		MinPalmScore = opt_min_palm_score;
	if (!opt_notrunclabels)
		opt_trunclabels = true;

	SetExcludes();

	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	const uint ThreadCount = GetRequestedThreadCount();
	vector<ObjMgr *> OMs;
	vector<RdRpSearcher *> Mods;
	for (int i = 0; i < int(ThreadCount); ++i)
		{
		ObjMgr *OM = new ObjMgr;
		OMs.push_back(OM);
		RdRpSearcher *RS = new RdRpSearcher;
		RS->Init(Model);
		Mods.push_back(RS);
		}

	FASTASeqSource *SS = new FASTASeqSource;
	SS->Open(QueryFileName);

	uint LastElapsedSecs = 0;
	uint CurrElapsedSecs = 0;
	ProgressStep(0, 1000, "Searching");
#pragma omp parallel num_threads(ThreadCount)
	{
	int ThreadIndex = GetThreadIndex();
	RdRpSearcher &RS = *Mods[ThreadIndex];
	RS.m_MinCScore = (float) MinCScore;
	RS.m_Trace = Trace;
	ObjMgr *OM = OMs[ThreadIndex];
	for (;;)
		{
		SeqInfo *QSI = OM->GetSeqInfo();

		if (ThreadIndex == 0)
			{
			if (g_QueryCount%100 == 0)
				{
				CurrElapsedSecs = GetElapsedSecs();
				if (CurrElapsedSecs > LastElapsedSecs)
					{
					uint Pct10 = SS->GetPctDoneX10();
					double HitPct = GetPct(g_FoundCount, g_QueryCount);
					ProgressStep(Pct10, 1000, "Searching %u / %u hits (%.1f%%)",
					  g_FoundCount, g_QueryCount, HitPct);
					LastElapsedSecs = CurrElapsedSecs;
					}
				}
			}

		bool Ok = SS->GetNext(QSI);
		if (!Ok)
			break;

		const string Label = string(QSI->m_Label);
		string Seq;
		for (uint i = 0; i < QSI->m_L; ++i)
			Seq += char(toupper(QSI->m_Seq[i]));

		RS.CFilter(Label, Seq);
		OM->Down(QSI);

		Lock();
		++g_QueryCount;
		Unlock();
		}
	}

	double HitPct = GetPct(g_FoundCount, g_QueryCount);
	ProgressStep(999, 1000, "Searching %u/%u hits (%.1f%%)",
	  g_FoundCount, g_QueryCount, HitPct);
	}

void cmd_searchc()
	{
	Die("Too many FPs");
	const string &QueryFileName = opt_searchc;
	SearchC(QueryFileName);
	}
