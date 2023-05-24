#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "outputfiles.h"

void cmd_pdb2cal()
	{
	const string &FN = opt_pdb2cal;

	ChainReader CR;
	CR.Open(FN);

	PDBChain Chain;
	uint Count = 0;
	for (;;)
		{
		bool Ok = CR.GetNext(Chain);
		if (!Ok)
			break;

		if (++Count%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u converted >%s\r",
			  sPct.c_str(), Count, Chain.m_Label.c_str());
			}
		Chain.ToCal(g_fcal);
		}
	Progress("100.0%% done, %u converted\n", Count);
	}

void cmd_pdb2cal_files()
	{
	const string &FN = opt_pdb2cal_files;

	vector<string> Paths;
	ReadLinesFromFile(FN, Paths);

	const uint N = SIZE(Paths);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Converting");

		const string &FileName = Paths[i];

		vector<PDBChain *> Chains;
		ReadChains(FileName, Chains);
		for (uint j = 0; j < SIZE(Chains); ++j)
			{
			Chains[j]->ToCal(g_fcal);
			delete Chains[j];
			}
		}
	}
