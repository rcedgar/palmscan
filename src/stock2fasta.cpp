#include "myutils.h"
#include "stock.h"
#include "outputfiles.h"

void cmd_stock2fasta()
	{
	const string &InputFN = opt_stock2fasta;

	Stock S;
	S.FromFile(InputFN);
	S.ToFasta(g_ffasta);
	}
