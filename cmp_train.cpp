#include "myutils.h"
#include "pdbchain.h"
#include "cmpsearcher.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "chainreader.h"
#include "cmp.h"
#include <map>

static uint g_DoneCount;
static uint g_HitCount;

static double g_MinScoreA = DBL_MAX;
static double g_MinScoreB = DBL_MAX;
static double g_MinScoreC = DBL_MAX;

static double g_MinScoreAB = DBL_MAX;
static double g_MinScoreBC = DBL_MAX;
static double g_MinScoreAC = DBL_MAX;

static double g_MinScore3 = DBL_MAX;

void CMP::GetMeanStdDev(uint i, uint j,
  double &Mean, double &StdDev) const
	{
	const uint N = SIZE(m_DistMxVec);
	asserta(N > 0);
	double Sum = 0;
	for (uint k = 0; k < N; ++k)
		{
		asserta(i < SIZE(m_DistMxVec[k]));
		asserta(j < SIZE(m_DistMxVec[k][i]));
		double d = m_DistMxVec[k][i][j];
		Sum += d;
		}
	Mean = double(Sum)/N;

	double Sumd2 = 0;
	for (uint k = 0; k < N; ++k)
		{
		double d = m_DistMxVec[k][i][j];
		double d2 = (d - Mean)*(d - Mean);
		Sumd2 += d2;
		}
	StdDev = (double) sqrt(Sumd2/N);
	}

void CMP::FinalizeTrain()
	{
	for (uint i = 0; i < PPSPL; ++i)
		{
		m_RefMeans[i][i] = 0;
		m_StdDevs[i][i] = 0;
		for (uint j = 0; j < i; ++j)
			{
			double Mean, StdDev;
			GetMeanStdDev(i, j, Mean, StdDev);

			m_RefMeans[i][j] = Mean;
			m_RefMeans[j][i] = Mean;

			m_StdDevs[i][j] = StdDev;
			m_StdDevs[j][i] = StdDev;
			}
		}
	}

void CMP::TrainChain(const PDBChain &Chain,
  uint APos, uint BPos, uint CPos)
	{
	asserta(APos != UINT_MAX);
	asserta(BPos != UINT_MAX);
	asserta(CPos != UINT_MAX);

	vector<vector<double> > DistMx;
	GetDistMx(Chain, APos, BPos, CPos, DistMx);

	m_RefLabels.push_back(Chain.m_Label);
	m_DistMxVec.push_back(DistMx);
	}

void CMP::ToDBFile(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);

	CloseStdioFile(f);
	}

static void Train(ChainReader &CR,
  const vector<string> &Labels, 
  const map<string, uint> &LabelToIndex, 
  const vector<vector<uint> > &MotifCoordsVec,
  CMP &Prof)
	{
	g_DoneCount = 0;
	g_HitCount = 0;

	g_MinScoreA = DBL_MAX;
	g_MinScoreB = DBL_MAX;
	g_MinScoreC = DBL_MAX;

	g_MinScoreAB = DBL_MAX;
	g_MinScoreBC = DBL_MAX;
	g_MinScoreAC = DBL_MAX;

	g_MinScore3 = DBL_MAX;

	PDBChain Q;
	uint NotFoundCount = 0;
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("Train %s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}

		const string &Seq = Q.m_Seq;
		const string Label = Q.m_Label;

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

		if (Seq[APos+3] != 'D')
			{
			Warning("Training error: MotifA in >%s\n", Label.c_str());
			Log("[%3u] ", APos+1);
			for (uint i = APos; i < APos + AL; ++i)
				Log("%c", Seq[i]);
			Log("\n");
			continue;
			}
		if (Seq[BPos+1] != 'G')
			{
			Warning("Training error: MotifB in >%s\n", Label.c_str());
			Log("[%3u] ", BPos+1);
			for (uint i = BPos; i < BPos + BL; ++i)
				Log("%c", Seq[i]);
			Log("\n");
			continue;
			}
		if (Seq[CPos+3] != 'D')
			{
			Warning("Training error: MotifC in >%s\n", Label.c_str());
			Log("[%3u] ", CPos+1);
			for (uint i = CPos; i < CPos + CL; ++i)
				Log("%c", Seq[i]);
			Log("\n");
			continue;
			}

		++g_HitCount;
		Prof.TrainChain(Q, APos, BPos, CPos);
		}
	Progress("Train 100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);
	if (NotFoundCount > 0)
		Warning("%u labels not found", NotFoundCount);
	asserta(g_HitCount > 0);
	Prof.FinalizeTrain();
	}

static void Write1(const char *Method, const PDBChain &Q,
  const CMPSearcher &CS, uint APos, uint BPos, uint CPos)
	{
	if (APos == UINT_MAX || BPos == UINT_MAX || CPos == UINT_MAX)
		{
		Log("%6.6s", ".");
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %12.12s", ".");
		Log("  %14.14s", ".");
		Log("  %8.8s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %6.6s", ".");
		Log("  %s", Method);
		Log("  >%s\n", Q.m_Label.c_str());
		return;
		}

	double ScoreA = CS.GetScoreA(Q, APos);
	double ScoreB = CS.GetScoreB(Q, BPos);
	double ScoreC = CS.GetScoreC(Q, CPos);

	double ScoreAB = CS.GetScoreAB(Q, APos, BPos);
	double ScoreBC = CS.GetScoreBC(Q, BPos, CPos);
	double ScoreAC = CS.GetScoreAC(Q, APos, CPos);

	g_MinScoreA = min(ScoreA, g_MinScoreA);
	g_MinScoreB = min(ScoreB, g_MinScoreB);
	g_MinScoreC = min(ScoreC, g_MinScoreC);

	g_MinScoreAB = min(ScoreAB, g_MinScoreAB);
	g_MinScoreBC = min(ScoreBC, g_MinScoreBC);
	g_MinScoreAC = min(ScoreAC, g_MinScoreAC);

	double Score3 = 0;

	Score3 += ScoreA;
	Score3 += ScoreB;
	Score3 += ScoreC;

	Score3 += ScoreAB;
	Score3 += ScoreBC;
	Score3 += ScoreAC;

	Score3 /= 6;

	g_MinScore3 = min(Score3, g_MinScore3);

	const char *Seq = Q.m_Seq.c_str();

	Log("%6.4f", Score3);
	Log("  %4u", APos + 1);
	Log("  %4u", BPos + 1);
	Log("  %4u", CPos + 1);
	Log("  %12.12s", Seq + APos);
	Log("  %14.14s", Seq + BPos);
	Log("  %8.8s", Seq + CPos);
	Log("  %6.4f", ScoreA);
	Log("  %6.4f", ScoreB);
	Log("  %6.4f", ScoreC);
	Log("  %6.4f", ScoreAB);
	Log("  %6.4f", ScoreBC);
	Log("  %6.4f", ScoreAC);
	Log("  %s", Method);
	Log("  >%s\n", Q.m_Label.c_str());
	}

static void Test(ChainReader &CR,
  const vector<string> &Labels, 
  const map<string, uint> &LabelToIndex, 
  const vector<vector<uint> > &MotifCoordsVec,
  CMPSearcher &CS)
	{
	g_DoneCount = 0;
	g_HitCount = 0;

	PDBChain Q;

	Log("\n");
	Log("%6.6s", "Score");
	Log("  %4.4s", "PosA");
	Log("  %4.4s", "PosB");
	Log("  %4.4s", "PosC");
	Log("  %12.12s", "SeqA");
	Log("  %14.14s", "SeqB");
	Log("  %8.8s", "SeqC");
	Log("  %6.6s", "A");
	Log("  %6.6s", "B");
	Log("  %6.6s", "C");
	Log("  %6.6s", "AB");
	Log("  %6.6s", "BC");
	Log("  %6.6s", "AC");
	Log("         Label\n");
	for (;;)
		{
		bool Ok = CR.GetNext(Q);
		if (!Ok)
			break;

		if (++g_DoneCount%100 == 0)
			{
			string sPct;
			CR.GetStrPctDone(sPct);
			Progress("Test %s%% done, %u / %u hits\r",
			  sPct.c_str(), g_HitCount, g_DoneCount);
			}

		const string &Seq = Q.m_Seq;
		string Label = Q.m_Label;

		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			continue;
		uint Index = p->second;

		asserta(Index < SIZE(MotifCoordsVec));
		const vector<uint> &MotifCoords = MotifCoordsVec[Index];
		asserta(SIZE(MotifCoords) == 3);

		uint APos_Train = MotifCoords[0];
		uint BPos_Train = MotifCoords[1];
		uint CPos_Train = MotifCoords[2];

		uint APos_PS = UINT_MAX;
		uint BPos_PS = UINT_MAX;
		uint CPos_PS = UINT_MAX;
		CS.Search(Q);
		CS.GetPSSMStarts(APos_PS, BPos_PS, CPos_PS);
		const char *Method = "???";
		if (APos_PS == APos_Train
			&& BPos_PS == BPos_Train
			&& CPos_PS == CPos_Train)
			Method = " Same";
		else
			Method = "*DIFF";
		Write1(Method, Q, CS, APos_PS, BPos_PS, CPos_PS);
		Write1("Train", Q, CS, APos_Train, BPos_Train, CPos_Train);

		++g_HitCount;
		Log("\n");
		}

	Log("\n");
	if (g_HitCount > 0)
		{
		Log("%6.4f", g_MinScore3);
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %4.4s", ".");
		Log("  %12.12s", ".");
		Log("  %14.14s", ".");
		Log("  %8.8s", ".");
		Log("  %6.4f", g_MinScoreA);
		Log("  %6.4f", g_MinScoreB);
		Log("  %6.4f", g_MinScoreC);
		Log("  %6.4f", g_MinScoreAB);
		Log("  %6.4f", g_MinScoreBC);
		Log("  %6.4f", g_MinScoreAC);
		Log("  Minimum\n");
		}

	Progress("Test 100.0%% done, %u / %u hits\r", g_HitCount, g_DoneCount);

	Log("%u / %u CMP hits\n", g_HitCount, g_DoneCount);
	}

void ReadMotifCoords(
  vector<vector<uint> > &MotifCoordsVec,
  vector<string> &Labels,
  map<string, uint> &LabelToIndex)
	{
	MotifCoordsVec.clear();
	Labels.clear();
	LabelToIndex.clear();

	if (!optset_motif_coords)
		Die("Must specify -motif_coords");

	FILE *f = OpenStdioFile(opt_motif_coords);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 7);

		const string &Label = Fields[0];

		uint LoA = StrToUint(Fields[1]);
		uint HiA = StrToUint(Fields[2]);
		asserta(LoA > 0);
		asserta(LoA < HiA);
		asserta(HiA - LoA + 1 == AL);

		uint LoB = StrToUint(Fields[3]);
		uint HiB = StrToUint(Fields[4]);
		asserta(LoB > 0);
		asserta(LoB < HiB);
		asserta(HiB - LoB + 1 == BL);

		uint LoC = StrToUint(Fields[5]);
		uint HiC = StrToUint(Fields[6]);
		asserta(LoC > 0);
		asserta(LoC < HiC);
		asserta(HiC - LoC + 1 == CL);

		vector<uint> MotifCoords;
		MotifCoords.push_back(LoA-1);
		MotifCoords.push_back(LoB-1);
		MotifCoords.push_back(LoC-1);

		uint Index = SIZE(Labels);
		asserta(LabelToIndex.find(Label) == LabelToIndex.end());
		LabelToIndex[Label] = Index;
		Labels.push_back(Label);
		MotifCoordsVec.push_back(MotifCoords);
		}
	}

void cmd_cmp_train()
	{
	const string &QueryFN = opt_cmp_train;

	vector<vector<uint> > MotifCoordsVec;
	vector<string> Labels;
	map<string, uint> LabelToIndex;
	ReadMotifCoords(MotifCoordsVec, Labels, LabelToIndex);

	PDBChain Q;

	ChainReader CR;
	CR.Open(QueryFN);

	CMP Prof;
	Train(CR, Labels, LabelToIndex, MotifCoordsVec, Prof);
	Prof.ToFile(opt_output);

	CR.Clear();
	CR.Open(QueryFN);

	CMPSearcher CS;
	CR.Clear();
	CR.Open(QueryFN);
	CS.m_DistMx = &Prof.m_RefMeans;
	CS.m_StdDevs = &Prof.m_StdDevs;
	Test(CR, Labels, LabelToIndex, MotifCoordsVec, CS);
	}
