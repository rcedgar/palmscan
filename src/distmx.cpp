#include "myutils.h"
#include "pdbchain.h"
#include "outputfiles.h"

void cmd_distmx()
	{
	const string &QFN = opt_distmx;
	asserta(g_ftsv != 0);

	vector<PDBChain *> Qs;
	ReadChains(QFN, Qs);

	const uint N = SIZE(Qs);
	
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Q = *Qs[i];
		uint QL = Q.GetSeqLength();
		const char *LabelQ = Q.m_Label.c_str();
		const string &SeqQ = Q.m_Seq;
		for (uint i = 0; i < QL; ++i)
			{
			char ci = SeqQ[i];
			for (uint j = i+1; j < QL; ++j)
				{
				char cj = SeqQ[j];
				double d = Q.GetDist(i, j);
				fprintf(g_ftsv, "%s", LabelQ);
				fprintf(g_ftsv, "\t%u", i);
				fprintf(g_ftsv, "\t%u", j);
				fprintf(g_ftsv, "\t%c", ci);
				fprintf(g_ftsv, "\t%c", cj);
				fprintf(g_ftsv, "\t%.1f", d);
				fprintf(g_ftsv, "\n");
				}
			}
		}
	}
