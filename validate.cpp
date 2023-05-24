#include "myutils.h"
#include "pdbchain.h"
#include "ppcsearcher.h"
#include "calreader.h"
#include "tshitmgr.h"
#include "outputfiles.h"

void SplitTrainTest(const vector<PDBChain *> &Chains, double TrainPct,
  vector<PDBChain *> &TrainChains, vector<PDBChain *> &TestChains);

static bool ChainToPPC(PDBChain &Q,
  const vector<PDBChain *> Refs, PDBChain &PPC)
	{
	PPC.Clear();

	TSHitMgr HM;
	HM.SetQuery(Q);

	const uint RefN = SIZE(Refs);
	TriSearcher TS;
	for (uint iR = 0; iR < RefN; ++iR)
		{
		const PDBChain &R = *Refs[iR];

		if (opt_self && Q.m_Label == R.m_Label)
			continue;

		TS.Search(Q, R);

		TSHit TH;
		bool Ok = TS.GetTopHit(TH);
		if (Ok)
			HM.Add(TH);
		}
	if (HM.m_Hits.empty())
		return false;

	HM.SetTopHit();
	HM.GetPPC(PPC);
	return true;
	}

static void Validate1(const vector<PDBChain *> &InputChains,
  double TrainPct, uint RefN, uint &TP, uint &FP, uint &TN, uint &FN)
	{
	TP = 0;
	FP = 0;
	TN = 0;
	FN = 0;

	vector<PDBChain *> TrainChains;
	vector<PDBChain *> TestChains;
	SplitTrainTest(InputChains, TrainPct, TrainChains, TestChains);

	const uint NTrain = SIZE(TrainChains);
	const uint NTest = SIZE(TestChains);

	if (RefN > NTrain)
		RefN = NTrain;

	vector<PDBChain *> RefChains;
	for (uint i = 0; i < RefN; ++i)
		{
		PDBChain *Chain = TrainChains[i];
		RefChains.push_back(Chain);
		}

	vector<PDBChain *> TestPPCs;
	uint PPCFound = 0;
	uint PPCFailed = 0;
	for (uint i = 0; i < NTest; ++i)
		{
		ProgressStep(i, NTest, "Convert to PPC");
		PDBChain &Chain = *TestChains[i];
		PDBChain *TestPPC = new PDBChain;
		bool Ok = ChainToPPC(Chain, RefChains, *TestPPC);
		if (Ok)
			{
			TestPPCs.push_back(TestPPC);
			++PPCFound;
			}
		else
			{
			++PPCFailed;
			TestPPCs.push_back(0);
			}
		}
	ProgressLog("%u PPC found, %u failed\n", PPCFound, PPCFailed);

	vector<bool> TrainIsRdRpVec;
	for (uint i = 0; i < NTrain; ++i)
		{
		const string &Label = TrainChains[i]->m_Label;
		bool IsRdRp = StartsWith(Label, "rdrp.");
		TrainIsRdRpVec.push_back(IsRdRp);
		}

	vector<bool> TestIsRdRpVec;
	for (uint i = 0; i < NTest; ++i)
		{
		const string &Label = TestChains[i]->m_Label;
		bool IsRdRp = StartsWith(Label, "rdrp.");
		TestIsRdRpVec.push_back(IsRdRp);
		}

	PPCSearcher PS;
	PS.SetRefs(RefChains);
	PS.InitClassify();

	uint RecCount = 0;
	uint HitCount = 0;
	PDBChain Q;
	for (uint TestIndex = 0; TestIndex < NTest; ++TestIndex)
		{
		bool TestIsRdRp = TestIsRdRpVec[TestIndex];
		PDBChain *ptrPPC = TestPPCs[TestIndex];
		if (ptrPPC == 0)
			{
			if (TestIsRdRp)
				++FN;
			else
				++TN;
			continue;
			}

		PDBChain &Q = *ptrPPC;
		PS.Search(Q);
		PS.WriteClassify(g_fclassify_tsv);
		asserta(PS.m_Classified);

		bool PredRdRp = PS.m_ClassifiedAsRdRp;
		if (TestIsRdRp && PredRdRp)
			++TP;
		else if (TestIsRdRp && !PredRdRp)
			++FN;
		else if (!TestIsRdRp && PredRdRp)
			++FP;
		else if (!TestIsRdRp && !PredRdRp)
			++TN;
		else
			asserta(false);
		}
	}

static void ValidateL1O(const vector<PDBChain *> &InputChains,
  uint RefN, uint &TP, uint &FP, uint &TN, uint &FN)
	{
	TP = 0;
	FP = 0;
	TN = 0;
	FN = 0;

	const uint N = SIZE(InputChains);
	if (RefN > N)
		RefN = N;

	vector<PDBChain *> RefChains;
	for (uint i = 0; i < RefN; ++i)
		{
		PDBChain *Chain = InputChains[i];
		RefChains.push_back(Chain);
		}

	vector<PDBChain *> TestPPCs;
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Convert to PPC");
		PDBChain &Chain = *InputChains[i];
		PDBChain *TestPPC = new PDBChain;
		bool Ok = ChainToPPC(Chain, RefChains, *TestPPC);
		if (Ok)
			TestPPCs.push_back(TestPPC);
		else
			TestPPCs.push_back(0);
		}

	vector<bool> TestIsRdRpVec;
	for (uint i = 0; i < N; ++i)
		{
		const string &Label = InputChains[i]->m_Label;
		bool IsRdRp = StartsWith(Label, "rdrp.");
		TestIsRdRpVec.push_back(IsRdRp);
		}

	PPCAligner PA;

	uint RecCount = 0;
	uint HitCount = 0;
	PDBChain Q;
	for (uint TestIndex = 0; TestIndex < N; ++TestIndex)
		{
		ProgressStep(TestIndex, N, "Leave-one-out");
		bool TestIsRdRp = TestIsRdRpVec[TestIndex];
		PDBChain *ptrPPC = TestPPCs[TestIndex];
		if (ptrPPC == 0)
			{
			if (TestIsRdRp)
				++FN;
			else
				++TN;
			continue;
			}

		PDBChain &Q = *ptrPPC;
		PA.SetQuery(Q);

		double TopRMSD = DBL_MAX;
		bool TopIsRdRp = false;
		for (uint TrainIndex = 0; TrainIndex < N; ++TrainIndex)
			{
			const PDBChain &R = *InputChains[TrainIndex];
			if (TrainIndex == TestIndex)
				continue;
			PA.SetRef(R);
			double RMSD = PA.GetMotifRMSD();
			if (RMSD < TopRMSD)
				{
				TopRMSD = RMSD;
				TopIsRdRp = TestIsRdRpVec[TrainIndex];
				}
			}
		if (TestIsRdRp && TopIsRdRp)
			++TP;
		else if (TestIsRdRp && !TopIsRdRp)
			++FN;
		else if (!TestIsRdRp && TopIsRdRp)
			++FP;
		else if (!TestIsRdRp && !TopIsRdRp)
			++TN;
		else
			asserta(false);
		}
	}

void cmd_validate()
	{
	const string &InputFN = opt_validate;

	asserta(optset_train_pct || optset_leave_one_out);
	asserta(optset_refn);

	uint RefN = opt_refn;
	double TrainPct = opt_train_pct;
	if (optset_leave_one_out)
		TrainPct = -1.0;

	vector<PDBChain *> InputChains;
	ReadChains(InputFN, InputChains);
	random_shuffle(InputChains.begin(), InputChains.end());

	const uint N = SIZE(InputChains);
	for (uint i = 0; i < N; ++i)
		{
		PDBChain &Q = *InputChains[i];
		Q.SetSS();
		}

	uint TP, FP, TN, FN;
	if (opt_leave_one_out)
		{
		ProgressLog("RefN %u, leave-one-out\n", RefN, TrainPct);
		ValidateL1O(InputChains, RefN, TP, FP, TN, FN);
		}
	else
		{
		ProgressLog("RefN %u, TrainPct %.1f%%\n", RefN, TrainPct);
		Validate1(InputChains, TrainPct, RefN, TP, FP, TN, FN);
		}
	ProgressLog("TP %u, FP %u, TN %u, FN %u\n", TP, FP, TN, FN);
	if (g_ftsv != 0)
		{
		fprintf(g_ftsv, "%u", RefN);
		if (opt_leave_one_out)
			fprintf(g_ftsv, "\tleave-one-out");
		else
			fprintf(g_ftsv, "\t%.1f", TrainPct);
		fprintf(g_ftsv, "\t%u\t%u\t%u\t%u", TP, FP, TN, FN);
		fprintf(g_ftsv, "\n");
		}
	}
