#include "myutils.h"
#include "pssm.h"
#include "rdrpsearcher.h"
#include "seqinfo.h"
#include "fastaseqsource.h"
#include "omplock.h"
#include <time.h>

static uint g_QueryCount;
static uint g_FoundCount;

vector<string> g_ExcludeNames;
vector<string> g_IncludeNames;

void SetExcludes()
	{
	if (!optset_exclude)
		return;

	string NamesStr = string(opt_exclude);
	Split(NamesStr, g_ExcludeNames, '+');
	const uint N = SIZE(g_ExcludeNames);
	for (uint i = 0; i < N; ++i)
		ProgressLog("  Exclude %s\n", g_ExcludeNames[i].c_str());
	}

void SetIncludes()
	{
	if (!optset_include)
		return;

	string NamesStr = string(opt_include);
	Split(NamesStr, g_IncludeNames, '+');
	const uint N = SIZE(g_IncludeNames);
	for (uint i = 0; i < N; ++i)
		ProgressLog("  Include %s\n", g_IncludeNames[i].c_str());
	}

void SeqToUpper(string &Seq)
	{
	const uint QL = SIZE(Seq);
	for (uint i = 0; i < QL; ++i)
		Seq[i] = toupper(Seq[i]);
	}

void SearchPSSMs(const string &QueryFileName, fn_OnPalmHit OnHit)
	{
	bool Trace = opt_trace;
	double MinPalmScore = 10.0;
	if (optset_min_palm_score)
		MinPalmScore = opt_min_palm_score;

	if (!opt_notrunclabels)
		opt_trunclabels = true;

	SetExcludes();

	RdRpModel Model;
	if (optset_model)
		{
		const string &ModelFileName = opt_model;
		Model.FromModelFile(ModelFileName);
		}
	else
		{
		vector<string> Lines;
		Model.GetDefaultModelLines(Lines);
		Model.FromLines(Lines);
		}

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

	//RdRpSearcher::InitOutput();
	uint LastElapsedSecs = 0;
	uint CurrElapsedSecs = 0;
	ProgressStep(0, 1000, "Searching");
#pragma omp parallel num_threads(ThreadCount)
	{
	int ThreadIndex = GetThreadIndex();
	RdRpSearcher &RS = *Mods[ThreadIndex];
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

		RS.Search(Label, Seq);
		bool IsHit = (RS.m_TopPalmHit.m_Score >= MinPalmScore);
		if (IsHit)
			{
			RS.WriteOutput();
#pragma omp critical
			{
			if (OnHit != 0)
				(*OnHit)(RS);
			++g_FoundCount;
			}
			}
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

void cmd_search_pssms()
	{
	const string &QueryFileName = opt_search_pssms;
	SearchPSSMs(QueryFileName, 0);
	}
