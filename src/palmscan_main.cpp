#include "myutils.h"
#include "timing.h"

int g_Frame = 0;

int main(int argc, char **argv)
	{
	MyCmdLine(argc, argv);
	LogProgramInfoAndCmdLine();
	if (!opt_quiet)
		{
		PrintProgramInfo(stdout);
		PrintCopyright(stdout);
		}

	if (optset_frame)
		{
		g_Frame = atoi(opt_frame);
		asserta(g_Frame >= -3 && g_Frame <= 3 && g_Frame != 0);
		}

	if (0)
		;
#define x(opt, f)	else if (optset_##opt) { void f(); f(); }
	x(build_pssm, BuildPSM)
	x(search_aa_top, SearchAATop)
	x(search_nt_top, SearchNtTop)
	x(search_nt_topx, SearchNtTopX)
	x(search_nt_all, SearchNtAll)
	x(test, Test_Viterbi_PSSM)
	x(search_pp, SearchPP)
	x(build_rdrp_model, BuildRdRpModel)
	x(otumaps, OTUMaps)
	x(alnqc, AlnQC)
#undef x
	else
		Die("No command specified");

	//LogTiming();
	//LogAllocs();

	LogElapsedTimeAndRAM();
	return 0;
	}
