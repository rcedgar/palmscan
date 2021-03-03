#include "myutils.h"
#include "pssm.h"

void BuildPSM()
	{
	const string InputFileName = string(opt_build_pssm);
	const string OutputFileName = string(opt_output);
	if (OutputFileName == "")
		Die("Missing -output filename");

	PSSM P;
	P.FromFasta(InputFileName);
//	P.LogMe();
	P.ToFile(OutputFileName);
	
	PSSM P2;
	P2.FromFile(OutputFileName);
	P2.LogMe();
	}
