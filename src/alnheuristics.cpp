#include "myutils.h"
#include "alnparams.h"
#include "alnheuristics.h"
#include "alpha.h"

AlnHeuristics::AlnHeuristics()
	{
	BandRadius = 0;
	HSPFinderWordLength = 0;
	XDropG = 0.0f;
	XDropU = 0.0f;
	MinGlobalHSPLength = 0;
	}

void AlnHeuristics::LogMe() const
	{
	Log("AH: Band %u, HSPw %u, Xg %.1f, Xu %.1f, HSP %u\n",
	  BandRadius,
	  HSPFinderWordLength,
	  XDropG,
	  XDropU,
	  MinGlobalHSPLength);
	}

void AlnHeuristics::InitFromCmdLine(const AlnParams &AP)
	{
	XDropU = (float) opt_xdrop_u;
	XDropG = (float) opt_xdrop_g;
	XDropGlobalHSP = (float) opt_xdrop_nw;

	BandRadius = opt_band;
	MinGlobalHSPLength = opt_minhsp;

	if (AP.GetIsNucleo())
		{
		HSPFinderWordLength = 5;
		MinGlobalHSPFractId = 0.75f;;
		MinGlobalHSPScore = MinGlobalHSPFractId*MinGlobalHSPLength*(float) opt_match;
		}
	else
		{
		HSPFinderWordLength = 3;
		
	// Avg BLOSUM62 score on the diagonal is 5.2, for comparison
		const float * const *SubstMx = AP.SubstMx;
		float MinDiagScore = 9e9f;
		for (unsigned i = 0; i < 20; ++i)
			{
			byte c = g_LetterToCharAmino[i];
			float Score = SubstMx[c][c];
			if (Score < MinDiagScore)
				MinDiagScore = Score;
			}

		MinGlobalHSPFractId = 0.5f;
		MinGlobalHSPScore = MinGlobalHSPFractId*MinDiagScore*MinGlobalHSPLength;
		}

	if (optset_hspw)
		HSPFinderWordLength = opt_hspw;

	if (opt_fulldp)
		{
		InitGlobalFullDP();
		return;
		}
	}

void AlnHeuristics::InitGlobalFullDP()
	{
	MinGlobalHSPLength = 0;
	HSPFinderWordLength = 0;
	BandRadius = 0;
	}
