#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "trifinder.h"

void cmd_find_tri()
	{
	const string &FN = opt_find_tri;

	ChainReader CR;
	CR.Open(FN);

	TriFinder TF;
	TF.SetSigmas(2.0);

	bool DGD = opt_dgd;

	PDBChain Chain;
	uint Count = 0;
	for (;;)
		{
		bool Ok = CR.GetNext(Chain);
		uint Mil = CR.GetMilDone();
		ProgressStep(Mil, 1001, "Working");
		if (!Ok)
			break;

		TF.Find(Chain, DGD);
		TF.LogMe();
		}
	}
