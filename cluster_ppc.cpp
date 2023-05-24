#include "myutils.h"
#include "pdbchain.h"
#include "ppcaligner.h"
#include "outputfiles.h"
#include <list>

void ReadPpc(const string &FN, vector<PDBChain *> &Chains);

void cmd_cluster_ppc()
	{
	const string &InputFN = opt_cluster_ppc;

	asserta(optset_rmsd);
	const double MaxRMSD = opt_rmsd;

	vector<PDBChain *> Chains;
	ReadPpc(InputFN, Chains);

	const uint N = SIZE(Chains);

	PPCAligner PA;

	list<uint> Pending;
	for (uint i = 0; i < N; ++i)
		Pending.push_back(i);

	uint DoneCount = 0;
	uint ClusterCount = 0;
	uint MemberCount = 0;
	for (;;)
		{
		++ClusterCount;
		ProgressStep(DoneCount, N+2, "%u clusters, %u members",
		  ClusterCount, MemberCount);
		if (Pending.size() == 0)
			break;
		uint CentroidIndex = Pending.front();
		++DoneCount;
		asserta(CentroidIndex < SIZE(Chains));
		Pending.pop_front();

		const PDBChain &CentroidChain = *Chains[CentroidIndex];
		const string &CentroidLabel = CentroidChain.m_Label;
		PA.SetRef(CentroidChain);
		if (g_fppc != 0)
			CentroidChain.ToCal(g_fppc);

		vector<list<uint>::const_iterator> MemberIndexes;
		for (list<uint>::const_iterator p = Pending.begin();
		  p != Pending.end(); ++p)
			{
			uint Index = *p;
			asserta(Index < SIZE(Chains));
			PDBChain &Chain = *Chains[Index];
			PA.SetQuery(Chain);
			double RMSD = PA.GetMotifRMSD();
			if (RMSD <= MaxRMSD)
				{
				++MemberCount;
				if (g_ftsv != 0)
					{
					const string &MemberLabel = Chain.m_Label;
					fprintf(g_ftsv, "%s", CentroidLabel.c_str());
					fprintf(g_ftsv, "\t%s", MemberLabel.c_str());
					fprintf(g_ftsv, "\t%.2f", RMSD);
					fprintf(g_ftsv, "\n");
					}
				MemberIndexes.push_back(p);
				++DoneCount;
				}
			}

		for (vector<list<uint>::const_iterator>::const_iterator q = MemberIndexes.begin();
		  q != MemberIndexes.end(); ++q)
			Pending.erase(*q);
		}
	ProgressStep(N+1, N+2, "%u clusters, %u members",
		ClusterCount, MemberCount);
	Log("%.3f RMSD, %u clusters, %u members",
	  MaxRMSD, ClusterCount, MemberCount);
	}
