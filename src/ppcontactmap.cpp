#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "cmp.h"

static uint g_DoneCount;
static uint g_HitCount;

static uint GetSeqPos(uint i, uint APos, uint BPos, uint CPos)
	{
	if (i >= AL + BL)
		return CPos + i - (AL + BL);
	if (i >= AL)
		return BPos + i - AL;
	return APos + i;
	}

void PPContactMap(const PDBChain &Q, 
  uint APos, uint BPos, uint CPos)
	{
	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
		return;

	++g_HitCount;
	const uint N = AL + BL + CL;
	vector<vector<double> > Mx;
	Mx.resize(N);
	for (uint i = 0; i < N; ++i)
		Mx[i].resize(N, DBL_MAX);

	for (uint i = 0; i < N; ++i)
		{
		uint SeqPosi = GetSeqPos(i, APos, BPos, CPos);
		for (uint j = 0; j < N; ++j)
			{
			uint SeqPosj = GetSeqPos(j, APos, BPos, CPos);
			double d = Q.GetDist(SeqPosi, SeqPosj);
			Mx[i][j] = d;
			}
		}

	FILE *f = g_fppcontactmap_tsv;
	if (f != 0)
		{
		const char *Label = Q.m_Label.c_str();
#pragma omp critical
		for (uint i = 1; i < N; ++i)
			{
			fprintf(f, "%s\t%u", Label, i);
			for (uint j = 0; j <= i; ++j)
				fprintf(f, "\t%.6g", Mx[i][j]);
			fprintf(f, "\n");
			}
		}
	}

static void Thread(ChainReader &CR,
  const vector<vector<uint> > &MotifCoordsVec,
  const vector<string> &Labels,
  const map<string, uint> &LabelToIndex)
	{
	uint ThreadIndex = GetThreadIndex();
	PDBChain Q;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			return;
#pragma omp critical
		{
		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("%s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}
		}

		const string &Seq = Q.m_Seq;
		const string &Label = Q.m_Label;

		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			{
			Log("Label not found >%s\n", Label.c_str());
			continue;
			}
		uint Index = p->second;
		asserta(Index < SIZE(MotifCoordsVec));
		const vector<uint> &MotifCoords = MotifCoordsVec[Index];
		asserta(SIZE(MotifCoords) == 3);

		uint APos = MotifCoords[0];
		uint BPos = MotifCoords[1];
		uint CPos = MotifCoords[2];

		PPContactMap(Q, APos, BPos, CPos);
		}
	}

void cmd_ppcontactmap()
	{
	const string &QueryFN = opt_ppcontactmap;

	vector<vector<uint> > MotifCoordsVec;
	vector<string> Labels;
	map<string, uint> LabelToIndex;
	ReadMotifCoords(MotifCoordsVec, Labels, LabelToIndex);

	uint ThreadCount = GetRequestedThreadCount();

	vector<vector<PDBChain *> > QVecs(ThreadCount);
	
	uint ThreadFinishedCount = 0;
	uint DoneCount = 0;
	uint HitCount = 0;

	ChainReader CR;
	CR.Open(QueryFN);

#pragma omp parallel num_threads(ThreadCount)
	Thread(CR, MotifCoordsVec, Labels, LabelToIndex);

	Progress("100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	}
