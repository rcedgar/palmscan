#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"

void cmd_checkppc()
	{
	const string &InputFN = opt_checkppc;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Chain = *Chains[i];
		asserta(Chain.m_MotifPosVec.size() == 3);
		asserta(Chain.CheckPPCMotifCoords());
		}
	ProgressLog("%u PPC records checked ok\n", N);
	}
