#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"

void cmd_pdbinfo()
	{
	const string &InputFileName = opt_pdbinfo;

	ChainReader CR;
	CR.Open(InputFileName);

	PDBChain Chain;
	uint ChainCount = 0;
	while (CR.GetNext(Chain))
		{
		uint SeqLength = Chain.GetSeqLength();
		ProgressLog("%5u atoms >%s\n",
		  SeqLength,
		  Chain.m_Label.c_str());
		++ChainCount;
		}
	ProgressLog("%u chains\n", ChainCount);
	}
