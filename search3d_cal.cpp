#include "myutils.h"
#include "pdbchain.h"
#include "trisearcher.h"
#include "abcxyz.h"
#include "sort.h"
#include "tshitmgr.h"
#include "calreader.h"
#include "omplock.h"

void ReadPpc(const string &FN, vector<PDBChain *> &Chains);
void Search1(TriSearcher &TS, TSHitMgr &HM,
  PDBChain &Q, vector<PDBChain *> &RefPDBs);

static uint g_DoneCount;
static uint g_HitCount;

static void Thread(CalReader &CR, vector<PDBChain> &Qs,
  vector<PDBChain *> &RefPDBs, vector<TriSearcher> &TSs,
  vector<TSHitMgr> &HMs)
	{
	uint ThreadIndex = GetThreadIndex();
	PDBChain& Q = Qs[ThreadIndex];
	TriSearcher &TS = TSs[ThreadIndex];
	TSHitMgr &HM = HMs[ThreadIndex];

	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			return;
		Lock("Done");
		++g_DoneCount;
		Unlock("Done");
		if (g_DoneCount%100 == 0)
			{
			Lock("Progress");
			string sPct, sPct2;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits (%s)\r",
			  sPct.c_str(), g_HitCount, g_DoneCount, GetPctStr(g_HitCount, g_DoneCount, sPct2));
			Unlock("Progress");
			}
		TriSearcher &TS = TSs[ThreadIndex];
		TSHitMgr &HM = HMs[ThreadIndex];

		Q.SetSS();
		Search1(TS, HM, Q, RefPDBs);
		if (!HM.m_Hits.empty())
			++g_HitCount;
		}
	}

void cmd_search3d_cal()
	{
	const string &QueryFN = opt_search3d_cal;
	const string &RefFN = opt_ref;

	CalReader CR;
	CR.Open(QueryFN);

	vector<PDBChain *> RefPDBs;
	ReadPpc(RefFN, RefPDBs);
	for (uint i = 0; i < SIZE(RefPDBs); ++i)
		{
		PDBChain &R = *RefPDBs[i];
		R.SetSS();
		}

	const uint RefN = SIZE(RefPDBs);

	uint ThreadCount = GetRequestedThreadCount();
	vector<vector<PDBChain *> > QVecs(ThreadCount);
	vector<TriSearcher> TSs(ThreadCount);
	vector<TSHitMgr> HMs(ThreadCount);
	vector<PDBChain> Qs(ThreadCount);

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Qs, RefPDBs, TSs, HMs);

	Progress("%u done, %u hits\n", g_DoneCount, g_HitCount);
	}
