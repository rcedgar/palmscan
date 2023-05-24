#include "myutils.h"
#include "pamerger.h"
#include "outputfiles.h"
#include "seqdb.h"
#include <time.h>
#include <map>

void cmd_pamerge()
	{
	const string &InputFileName = opt_pamerge;
	FILE *fIn = OpenStdioFile(InputFileName);

	PAMerger PAM;

	double MinScoreRdRp = 75;
	double MaxScoreNotRdRp = 25;

	if (optset_minscore)
		MinScoreRdRp = opt_minscore;

	string Line;
	uint RecCount = 0;
	uint64 FileSize = GetStdioFileSize64(fIn);
	uint N100 = 0;
	uint N90 = 0;
	uint N50 = 0;
	uint N1 = 0;
	uint N0 = 0;
	time_t t0 = time(0);
	while (ReadLineStdioFile(fIn, Line))
		{
		++RecCount;
		if (RecCount%10000 == 0)
			{
			time_t t1 = time(0);
			if (t1 != t0)
				{
				t0 = t1;
				uint64 FilePos = GetStdioFilePos64(fIn);
				double Pct = GetPct(double(FilePos), double(FileSize));
				Progress("%u recs, %.1f%% done\r", RecCount, Pct);
				}
			}
			  
		PAM.FromLine(Line);
		PAM.ToFev(g_ffev);
		PAM.ToTsv(g_ftsv);
		double RdRpScore = PAM.GetRdRpScore();
		if (RdRpScore >= MinScoreRdRp)
			PAM.ToFasta(g_ffasta);

		if (RdRpScore == 100)
			++N100;
		else if (RdRpScore >= 90)
			++N90;
		else if (RdRpScore >= 50)
			++N50;
		else if (RdRpScore >= 1)
			++N1;
		else
			{
			asserta(RdRpScore < 1 && RdRpScore >= 0);
			++N0;
			}
		}
	Progress("%u recs, 100.0%% done\n", RecCount);
	ProgressLog("%12u   Sequences\n", RecCount);
	ProgressLog("%12u   RdRp score == 100 (%.1f%%)\n",
	  N100, GetPct(N100, RecCount));
	ProgressLog("%12u   RdRp score >=  90 (%.1f%%)\n",
	  N90, GetPct(N90, RecCount));
	ProgressLog("%12u   RdRp score >=  50 (%.1f%%)\n",
	  N50, GetPct(N50, RecCount));
	ProgressLog("%12u   RdRp score >    0 (%.1f%%)\n",
	  N1, GetPct(N1, RecCount));
	ProgressLog("%12u   RdRp score =    0 (%.1f%%)\n",
	  N0, GetPct(N0, RecCount));
	}
