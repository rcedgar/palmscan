#include "myutils.h"
#include "pdbchain.h"
#include "chainreader.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include <map>

void ReadMotifCoords(
  vector<vector<uint> > &MotifCoordsVec,
  vector<string> &Labels,
  map<string, uint> &LabelToIndex);

void GetHSE(const PDBChain &Chain, uint Pos, double Radius,
  uint MinPos, uint MaxPos, uint &NU, uint &ND);

void GetHSE_NUs(const PDBChain &Chain, double Radius,
  uint APos, uint BPos, uint CPos,
  vector<uint> &NUs)
	{
	uint MinPos = APos;
	uint MaxPos = CPos + CL - 1;

	NUs.clear();

	for (uint Pos = APos; Pos < APos + AL; ++Pos)
		{
		uint NU, ND;
		GetHSE(Chain, Pos, Radius, MinPos, MaxPos, NU, ND);
		NUs.push_back(NU);
		}

	for (uint Pos = BPos; Pos < BPos + BL; ++Pos)
		{
		uint NU, ND;
		GetHSE(Chain, Pos, Radius, MinPos, MaxPos, NU, ND);
		NUs.push_back(NU);
		}

	for (uint Pos = CPos; Pos < CPos + CL; ++Pos)
		{
		uint NU, ND;
		GetHSE(Chain, Pos, Radius, MinPos, MaxPos, NU, ND);
		NUs.push_back(NU);
		}
	}

void cmd_hse_train()
	{
	const string &QueryFN = opt_hse_train;

	vector<vector<uint> > MotifCoordsVec;
	vector<string> Labels;
	map<string, uint> LabelToIndex;
	ReadMotifCoords(MotifCoordsVec, Labels, LabelToIndex);

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
		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			{
			if (NotFoundCount < 10)
				Log("Not found >%s\n", Label.c_str());
			++NotFoundCount;
			continue;
			}
		uint Index = p->second;

		asserta(Index < SIZE(MotifCoordsVec));
		const vector<uint> &MotifCoords = MotifCoordsVec[Index];
		asserta(SIZE(MotifCoords) == 3);

		uint APos = MotifCoords[0];
		uint BPos = MotifCoords[1];
		uint CPos = MotifCoords[2];

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
