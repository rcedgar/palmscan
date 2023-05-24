#include "myutils.h"
#include "pdbchain.h"
#include "trisearcher.h"
#include "abcxyz.h"
#include "sort.h"
#include "tshitmgr.h"
#include "omplock.h"

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);

void Search1(TriSearcher &TS, TSHitMgr &HM,
  PDBChain &Q, vector<PDBChain *> &RefPDBs)
	{
	HM.SetQuery(Q);

	const uint RefN = SIZE(RefPDBs);
	for (uint iR = 0; iR < RefN; ++iR)
		{
		const PDBChain &R = *RefPDBs[iR];

		if (opt_self && Q.m_Label == R.m_Label)
			continue;

		Q.SetSS();
		TS.Search(Q, R);
		TS.WriteOutput();
		TSHit TH;
		bool Ok = TS.GetTopHit(TH);
		if (Ok)
			HM.Add(TH);
		}

	HM.WriteOutput();
	}

void cmd_search3d_tri()
	{
	const string &QueryFN = opt_search3d_tri;
	const string &RefFN = opt_ref;

	vector<string> QueryFileNames;
	vector<string> RefFileNames;
	ReadLinesFromFile(QueryFN, QueryFileNames);
	ReadLinesFromFile(RefFN, RefFileNames);

	vector<PDBChain *> RefPDBs;
	ReadChains(RefFN, RefPDBs);
	for (uint i = 0; i < SIZE(RefPDBs); ++i)
		{
		PDBChain &Ref = *RefPDBs[i];
		if (Ref.m_MotifPosVec.size() != 3)
			Die("Reference missing motif spec: %s", Ref.m_Label.c_str());
		RefPDBs[i]->SetSS();
		}

	const uint QueryN = SIZE(QueryFileNames);
	const uint RefN = SIZE(RefPDBs);

	uint ThreadCount = GetRequestedThreadCount();
	vector<vector<PDBChain *> > QVecs(ThreadCount);
	vector<TriSearcher> TSs(ThreadCount);
	vector<TSHitMgr> HMs(ThreadCount);
	
	uint DoneCount = 0;

#pragma omp parallel for num_threads(ThreadCount)
	for (int iQ = 0; iQ < (int) QueryN; ++iQ)
		{
#pragma omp critical
		{
		ProgressStep(DoneCount++, QueryN, "Searching");
		}

		uint ThreadIndex = GetThreadIndex();
		vector<PDBChain *> &QVec = QVecs[ThreadIndex];
		TriSearcher &TS = TSs[ThreadIndex];
		TSHitMgr &HM = HMs[ThreadIndex];

		const string &QueryFileName = QueryFileNames[iQ];
		PDBChain::ReadChainsFromFile(QueryFileName, QVec);
		uint ChainCount = SIZE(QVec);
		for (uint ChainIndex = 0; ChainIndex < ChainCount; ++ChainIndex)
			{
			PDBChain &Q = *QVec[ChainIndex];
			Search1(TS, HM, Q, RefPDBs);
			}
		}
	}
