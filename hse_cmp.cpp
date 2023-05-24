#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "cmpsearcher.h"

void GetHSE(const PDBChain &Chain, uint Pos, double Radius,
  uint &NU, uint &ND);

void GetHSE_NUs(const PDBChain &Chain, double Radius,
  uint APos, uint BPos, uint CPos, vector<uint> &NUs);

void cmd_hse_cmp()
	{
	const string &QueryFN = opt_hse_cmp;

	if (!optset_model)
		Die("Must specify -model");
	const string &ModelFileName = opt_model;

	CMP Prof;
	Prof.FromFile(ModelFileName);

	CMPSearcher CS;
	CS.SetProf(Prof);

	PDBChain Chain;
	ChainReader CR;
	CR.Open(QueryFN);
	double Radius = 12.0;
	if (optset_radius)
		Radius = opt_radius;
	uint NotFoundCount = 0;
	vector<vector<uint> > NUVec;
	vector<string> FoundLabels;
	const uint ML = AL + BL + CL;
	for (;;)
		{
		bool Ok = CR.GetNext(Chain);
		if (!Ok)
			break;

		const string &Label = Chain.m_Label;
		CS.Search(Chain);
		uint APos, BPos, CPos;
		double PalmScore = CS.GetPSSMStarts(APos, BPos, CPos);
		if (PalmScore == 0)
			continue;

		Log("%u %u %u >%s\n", APos, BPos, CPos, Label.c_str());

		vector<uint> NUs;
		GetHSE_NUs(Chain, Radius, APos, BPos, CPos, NUs);
		asserta(SIZE(NUs) == ML);
		FoundLabels.push_back(Label);
		NUVec.push_back(NUs);
		}

	if (g_ftsv == 0)
		return;

	const uint N = SIZE(FoundLabels);
	asserta(SIZE(NUVec) == N);
	fprintf(g_ftsv, "Pos");
	for (uint i = 0; i < N; ++i)
		fprintf(g_ftsv, "\t%s", FoundLabels[i].c_str());
	fprintf(g_ftsv, "\n");

	for (uint Pos = 0; Pos < ML; ++Pos)
		{
		fprintf(g_ftsv, "%u", Pos);
		for (uint i = 0; i < N; ++i)
			fprintf(g_ftsv, "\t%u", NUVec[i][Pos]);
		fprintf(g_ftsv, "\n");
		}
	}
