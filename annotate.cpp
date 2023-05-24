#include "myutils.h"
#include "pdbchain.h"
#include "rdrpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "cmpsearcher.h"

void FindCavity(const PDBChain &Chain, const string &PDBFileName);

void cmd_annotate()
	{
	const string &QueryFN = opt_annotate;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	CMP Prof;
	Prof.FromFile(ModelFileName);

	ChainReader CR;
	CR.Open(QueryFN);

	CMPSearcher CS;
	CS.SetProf(Prof);

	PDBChain Q;
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

		PDBChain RotatedChain;
		Q.GetTriFormChain_DGD(RotatedChain);

		RotatedChain.ToPDB(opt_output);
		RotatedChain.ToPML(g_fpml, opt_output);

		FindCavity(RotatedChain, opt_output);
		break;
		}
	ProgressLog("%u chains converted\n", ConvertedCount);
	}
