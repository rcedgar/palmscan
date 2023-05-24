#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "quarts.h"
#include "outputfiles.h"
#include "trifinder.h"
#include "ppp.h"

void cmd_search_ppp()
	{
	const string &QueryFileName = opt_search_ppp;

	bool DGD = true;

	ChainReader CR;
	CR.Open(QueryFileName);

	uint ChainCount = 0;
	PDBChain Chain;
	PPP Profile;
	Profile.GlobalInit();
	Profile.FromPPP(opt_ppp);
	TriFinder TF;
	TF.SetSigmas(2);

	while (CR.GetNext(Chain))
		{
//		Chain.CheckPPCMotifCoords();
		uint SeqLength = Chain.GetSeqLength();
		uint Mil = CR.GetMilDone();
		ProgressStep(Mil, 1001, "Scanning");

		TF.Find(Chain, DGD);
		if (TF.m_Score <= 0)
			{
			Pf(g_ftsv, "0\t%s\n", Chain.m_Label.c_str());
			continue;
			}

		uint PosA = TF.m_PosA;
		uint PosB = TF.m_PosB;
		uint PosC = TF.m_PosC;
		double Score = Profile.ScoreChain(Chain, PosA, PosB, PosC);
		Pf(g_ftsv, "%.3g\t%s\n", Score, Chain.m_Label.c_str());
		}
	}
