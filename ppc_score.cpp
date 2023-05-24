#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "quarts.h"
#include "outputfiles.h"

double GetPPCProfileScore(
  uint LAB, uint LAC, uint LBC,
  double dAB, double dAC, double dBC);

void cmd_ppc_score()
	{
	const string &InputFileName = opt_ppc_score;

	ChainReader CR;
	CR.Open(InputFileName);

	Pf(g_ftsv, "Label");
	Pf(g_ftsv, "\tScore");
	Pf(g_ftsv, "\tDGD");
	Pf(g_ftsv, "\tPosA");
	Pf(g_ftsv, "\tPosB");
	Pf(g_ftsv, "\tPosC");
	Pf(g_ftsv, "\tLAB");
	Pf(g_ftsv, "\tLAC");
	Pf(g_ftsv, "\tLBC");
	Pf(g_ftsv, "\tdAB");
	Pf(g_ftsv, "\tdAC");
	Pf(g_ftsv, "\tdBC");
	Pf(g_ftsv, "\n");

	uint ChainCount = 0;
	PDBChain Chain;
	while (CR.GetNext(Chain))
		{
		uint SeqLength = Chain.GetSeqLength();
		if (++ChainCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% scanning profiles\r", sPct.c_str());
			}
		Chain.CheckPPCMotifCoords(true);

		uint PosA = Chain.m_MotifPosVec[A];
		uint PosB = Chain.m_MotifPosVec[B];
		uint PosC = Chain.m_MotifPosVec[C];

		char DGD[4];
		DGD[0] = Chain.m_Seq[PosA];
		DGD[1] = Chain.m_Seq[PosB];
		DGD[2] = Chain.m_Seq[PosC];
		DGD[3] = 0;

		double dAB = Chain.GetDist(PosA, PosB);
		double dAC = Chain.GetDist(PosA, PosC);
		double dBC = Chain.GetDist(PosB, PosC);

		uint LAB = PosB - PosA + 1;
		uint LAC = PosC - PosA + 1;
		uint LBC = PosC - PosB + 1;

		double Score = GetPPCProfileScore(LAB, LAC, LBC, dAB, dAC, dBC);

		Pf(g_ftsv, "%s", Chain.m_Label.c_str());
		Pf(g_ftsv, "\t%.3g", Score);
		Pf(g_ftsv, "\t%s", DGD);
		Pf(g_ftsv, "\t%u", PosA + 1);
		Pf(g_ftsv, "\t%u", PosB + 1);
		Pf(g_ftsv, "\t%u", PosC + 1);
		Pf(g_ftsv, "\t%u", LAB);
		Pf(g_ftsv, "\t%u", LAC);
		Pf(g_ftsv, "\t%u", LBC);
		Pf(g_ftsv, "\t%.2f", dAB);
		Pf(g_ftsv, "\t%.2f", dAC);
		Pf(g_ftsv, "\t%.2f", dBC);
		Pf(g_ftsv, "\n");
		}

	ProgressLog("%u chains\n", ChainCount);
	}
