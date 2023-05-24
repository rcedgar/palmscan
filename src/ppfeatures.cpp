#include "myutils.h"
#include "pdbchain.h"
#include "rdrpmodel.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"

#if 0

static uint g_DoneCount;
static uint g_HitCount;

static void LogSeqAndSS(const string &Msg, const string &Seq,
  const string SS, int Lo, int Hi, int CatPos)
	{
	Log("  %s[%3d]", Msg.c_str(), CatPos);

	int L = int(SIZE(Seq));

	if (Lo < 0)
		Lo = 0;
	if (Hi >= L)
		Hi = L - 1;

	string SeqRow;
	string SSRow;
	for (int i = Lo; i <= Hi; ++i)
		{
		if (i == CatPos)
			{
			SeqRow += ' ';
			SSRow += ' ';
			}

		SeqRow += Seq[i];
		SSRow += SS[i];

		if (i == CatPos)
			{
			SeqRow += ' ';
			SSRow += ' ';
			}
		}
	Log("  %s", SeqRow.c_str());
	Log("  %s", SSRow.c_str());
	}

static void PPFeatures(const PDBChain &Chain,
  uint APos, uint BPos, uint CPos)
	{
	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
		{
		//Log("Motif(s) not found in >%s\n", Chain.m_Label.c_str());
		return;
		}
	if (CPos < APos)
		{
		//Log("Permuted domain >%s\n", Chain.m_Label.c_str());
		return;
		}

	const string &Seq = Chain.m_Seq;

	string SS;
	Chain.GetSS(SS);

	vector<uint> PosVecA;
	vector<uint> PosVecB;
	vector<uint> PosVecC;
	SearchMotifA_D(Seq, SS, PosVecA);
	SearchMotifB_G(Seq, SS, PosVecB);
	SearchMotifC_D(Seq, SS, PosVecC);

	Log("\nAs: ");
	for (uint i = 0; i < SIZE(PosVecA); ++i)
		Log(" %u", PosVecA[i]);
	Log(", Bs: ");
	for (uint i = 0; i < SIZE(PosVecB); ++i)
		Log(" %u", PosVecB[i]);
	Log(", Cs: ");
	for (uint i = 0; i < SIZE(PosVecC); ++i)
		Log(" %u", PosVecC[i]);
	Log("\n");

	//Log("\n");
	//Log(">%s\n", Chain.m_Label.c_str());
	LogSeqAndSS("=A=", Seq, SS, int(APos) - 4, int(APos) + AL + 4, APos+3);
	LogSeqAndSS("=B=", Seq, SS, int(BPos) - 4, int(BPos) + CL + 4, BPos+1);
	LogSeqAndSS("=C=", Seq, SS, int(CPos) - 4, int(CPos) + AL + 4, CPos+3);
	Log("  >%s\n", Chain.m_Label.c_str());
	}

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

		PPFeatures(Q, APos, BPos, CPos);
		}
	}

void cmd_ppfeatures()
	{
	const string &QueryFN = opt_ppfeatures;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;
	RdRpModel Model;
	Model.FromModelFile(ModelFileName);

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
#endif
