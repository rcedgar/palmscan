#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"
#include <set>

void cmd_remove_dupes()
	{
	const string &InputFN = opt_remove_dupes;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);

	set<string> Seqs;
	uint DupeCount = 0;
	const uint N = SIZE(Chains);
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Chain = *Chains[i];
		const string &Seq = Chain.m_Seq;
		if (Seqs.find(Seq) != Seqs.end())
			{
			++DupeCount;
			continue;
			}
		Seqs.insert(Seq);
		Chain.ToCal(g_fcal);
		}
	ProgressLog("%u dupes deleted.\n", DupeCount);
	}
