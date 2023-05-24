#include "myutils.h"
#include "pdbchain.h"
#include "calreader.h"

void ReadPpc(const string &FN, vector<PDBChain *> &Chains)
	{
	Chains.clear();
	CalReader CR;
	CR.Open(FN);
	uint n = 0;
	for (;;)
		{
		++n;
		if (n%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("Reading PPCs %s (%s%%)\r", FN.c_str(), sPct.c_str());
			}
		PDBChain *Chain = new PDBChain;
		bool Ok = CR.GetNext(*Chain);
		if (!Ok)
			return;
		Chain->CheckPPCMotifCoords();
		Chains.push_back(Chain);
		}
	Progress("Reading PPCs %s (100.0%%)\n");
	}
