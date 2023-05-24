#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"

void cmd_cal2fa()
	{
	const string &InputFN = opt_cal2fa;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Chain = *Chains[i];
		const string &Seq = Chain.m_Seq;
		const string &Label = Chain.m_Label;
		if (!Seq.empty())
			SeqToFasta(g_ffasta, Label.c_str(), Seq.c_str(), SIZE(Seq));
		}
	}
