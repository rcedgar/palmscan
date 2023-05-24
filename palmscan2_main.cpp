#include "myutils.h"
#include "timing.h"
#include "outputfiles.h"

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

	extern vector<string> g_Argv;
	uint n = SIZE(g_Argv);
	asserta(n > 0);
	string ShortCmdLine = g_Argv[1];
	if (n > 2)
		ShortCmdLine += " " + g_Argv[2];

	ProgressPrefix(false);
	Progress("[%s]\n", ShortCmdLine.c_str() + 1);
	ProgressPrefix(true);

	OpenOutputFiles();
	if (0)
		;
#define C(x)	else if (optset_##x) { void cmd_##x(); cmd_##x(); }
#include "cmds.h"

	else
		Die("No command specified");

	CloseOutputFiles();

	LogElapsedTimeAndRAM();
	return 0;
	}
