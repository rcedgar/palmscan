#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "cmpsearcher.h"

void cmd_palmcore_cmp()
	{
	const string &QueryFN = opt_palmcore_cmp;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	CMP Prof;
	Prof.FromFile(ModelFileName);

	ChainReader CR;
	CR.Open(QueryFN);
	CMPSearcher CS;
	PDBChain Q;
	CS.SetProf(Prof);
	uint ConvertedCount = 0;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		CS.Search(Q);

		uint APos = UINT_MAX;
		uint BPos = UINT_MAX;
		uint CPos = UINT_MAX;
		double PalmScore = CS.GetPSSMStarts(APos, BPos, CPos);
		if (PalmScore <= 0)
			continue;

		++ConvertedCount;
		Q.SetMotifPosVec(APos, BPos, CPos);

		PDBChain PC;
		Q.GetPC(PC);

		PC.ToPDB(opt_output);
		PC.ToPML(g_fpml, opt_output);
		break;
		}
	ProgressLog("%u chains converted\n", ConvertedCount);
	}
