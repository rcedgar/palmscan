#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"

static uint g_DoneCount;
static uint g_HitCount;

static void Thread(ChainReader &CR, vector<PDBChain> &Qs,  
  vector<RdRpSearcher> &RSs)
	{
	uint ThreadIndex = GetThreadIndex();
	PDBChain &Q = Qs[ThreadIndex];
	RdRpSearcher &RS = RSs[ThreadIndex];
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			return;
#pragma omp critical
		{
		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}
		}

		const string &QSeq = Q.m_Seq;
		const string &QLabel = Q.m_Label;

		RS.Search(QLabel, QSeq);
		RS.WriteOutput();

		uint APos = RS.GetMotifPos(0);
		uint BPos = RS.GetMotifPos(1);
		uint CPos = RS.GetMotifPos(2);

		if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
			{
			Log("Motif(s) not found in >%s\n", QLabel.c_str());
			continue;
			}
		if (CPos < APos)
			{
			Log("Permuted domain >%s\n", QLabel.c_str());
			continue;
			}

#pragma omp critical
		{
		++g_HitCount;
		}

//		Q.m_MotifPosVec.clear();
//		Q.m_MotifPosVec.push_back(APos);
//		Q.m_MotifPosVec.push_back(BPos);
//		Q.m_MotifPosVec.push_back(CPos);
//
//		vector<vector<double> > MotifCoords;
//		Q.GetMotifCoords(MotifCoords);
//
//		vector<double> t;
//		vector<vector<double> > R;
//		GetTriForm(MotifCoords, t, R);
//
//		const uint QL = Q.GetSeqLength();
//		vector<double> Pt;
//		vector<double> XPt;
//		for (uint Pos = 0; Pos < QL; ++Pos)
//			{
//			Q.GetPt(Pos, Pt);
//			XFormPt(Pt, t, R, XPt);
//			Q.SetPt(Pos, XPt);
//			}
//		uint PPL = CPos + CL - APos;
//
//#pragma omp critical
//		{
//		PDBChain Q_PPC;
//		Q.SetMotifPosVec(APos, BPos, CPos);
//		Q.GetPPC(Q_PPC);
//		Q_PPC.ToCal(g_fppc);
//		if (optset_pdbout)
//			Q_PPC.ToPDB(opt_pdbout);
//		}
		}
	}

void cmd_search3d_pssms()
	{
	const string &QueryFN = opt_search3d_pssms;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

	//RdRpSearcher::InitOutput();
	uint ThreadCount = GetRequestedThreadCount();

	vector<vector<PDBChain *> > QVecs(ThreadCount);
	
	uint ThreadFinishedCount = 0;
	uint DoneCount = 0;
	uint HitCount = 0;
	vector<PDBChain> Qs(ThreadCount);
	vector<bool> ThreadDone(ThreadCount);
	vector<RdRpSearcher> RSs(ThreadCount);
	for (uint i = 0; i < ThreadCount; ++i)
		RSs[i].Init(Model);

	ChainReader CR;
	CR.Open(QueryFN);

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, Qs, RSs);

	Progress("100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	}
