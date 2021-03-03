#include "myutils.h"
#include "pssm.h"
#include "seqdb.h"
#include "rdrpmodel.h"

void Test()
	{
	const string &SpecFileName = opt_test;

	RdRpModel Mod;
	Mod.FromSpecFile(SpecFileName);
	}
