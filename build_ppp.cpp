#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "quarts.h"
#include "outputfiles.h"
#include "ppp.h"

void cmd_build_ppp()
	{
	const string &InputFileName = opt_build_ppp;

	FILE *fOut = CreateStdioFile(opt_output);

	ChainReader CR;
	CR.Open(InputFileName);

	uint ChainCount = 0;
	PDBChain Chain;
	PPP Profile;
	Profile.GlobalInit();
	Profile.Init(InputFileName);

	while (CR.GetNext(Chain))
		{
		uint SeqLength = Chain.GetSeqLength();
		uint Mil = CR.GetMilDone();
		ProgressStep(Mil, 1001, "Scanning");
		Profile.AddChain(Chain);
		}

	Profile.Averages();
	Profile.ToPPP(fOut);
	Profile.ToPDB(g_fpdb);
	CloseStdioFile(fOut);
	}
