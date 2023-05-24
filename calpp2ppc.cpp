#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"

// Input is CAL which is palmprint-only with coordinates,
//  output is rotated into PPC.
void cmd_calpp2ppc()
	{
	const string &InputFN = opt_calpp2ppc;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		ProgressStep(i, N, "Converting");
		const PDBChain &Chain = *Chains[i];
		Chain.CheckMotifCoords();

		PDBChain PPC;
		Chain.GetPPC(PPC);
		asserta(PPC.CheckPPCMotifCoords());
		PPC.ToCal(g_fppc);
		}
	}
