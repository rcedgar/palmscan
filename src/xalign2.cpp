#include "myutils.h"
#include "mx.h"
#include "xprof.h"
#include "xprofdata.h"
#include "xtrainer.h"
#include "outputfiles.h"
#include "omplock.h"

double GetDALIZ_PosVecs(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs);

void GetPosVecs(const string &QRow, const string &TRow,
  vector<uint> &PosQs, vector<uint> &PosTs,
  vector<uint> &MatchColToCol, vector<uint> &ColToMatchCol);

double GetDALIScore_LSB_Pair(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  uint Lo1, uint Hi1, uint Lo2, uint Hi2);

double GetDALIScore_Cols_Band(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  uint Radius, vector<double> &ColScores);

double GetDALIScore_LSB_Pair(const PDBChain &Q, const PDBChain &T,
  const vector<uint> &PosQs, const vector<uint> &PosTs,
  uint Lo1, uint Hi1, uint Lo2, uint Hi2);

static void GetLSBs(const vector<double> &Scores,
  uint MinLSBLength, vector<uint> &Los, vector<uint> &His)
	{
	Los.clear();
	His.clear();

	const uint L = SIZE(Scores);
	uint LSBLength = 0;
	uint Lo = UINT_MAX;
	for (uint i = 0; i < L; ++i)
		{
		double Score = Scores[i];
		if (Score >= 0)
			{
			if (Lo == UINT_MAX)
				Lo = i;
			}
		else
			{
			if (Lo != UINT_MAX && i - Lo >= MinLSBLength)
				{
				Los.push_back(Lo);
				His.push_back(i-1);
				}
			Lo = UINT_MAX;
			}
		}

	if (Lo != UINT_MAX && L - Lo >= MinLSBLength)
		{
		Los.push_back(Lo);
		His.push_back(L-1);
		}
	}

void XAlign2(
  const XProfData &ProfQ, const XProfData &ProfT,
  const string &RowQ, const string &RowT)
	{
	vector<uint> PosQs;
	vector<uint> PosTs;
	vector<uint> MatchColToCol;
	vector<uint> ColToMatchCol;
	GetPosVecs(RowQ, RowT, PosQs, PosTs, MatchColToCol, ColToMatchCol);

	vector<double> ColScores;
	GetDALIScore_Cols_Band(ProfQ, ProfT, PosQs, PosTs, 48, ColScores);

	vector<uint> PvePosQs;
	vector<uint> PvePosTs;
	const uint MatchColCount = SIZE(ColScores);
	asserta(SIZE(PosQs) == MatchColCount);
	asserta(SIZE(PosTs) == MatchColCount);
	for (uint i = 0; i < MatchColCount; ++i)
		{
		if (ColScores[i] > 0)
			{
			PvePosQs.push_back(PosQs[i]);
			PvePosTs.push_back(PosTs[i]);
			}
		}
	double Z = GetDALIZ_PosVecs(ProfQ, ProfT, PvePosQs, PvePosTs);
	Log("Z = %.3g\n", Z);
	return;

	vector<uint> Los;
	vector<uint> His;
	GetLSBs(ColScores, 16, Los, His);

	const uint LSBCount = SIZE(Los);
	asserta(SIZE(His) == LSBCount);
	for (uint i = 0; i < LSBCount; ++i)
		{
		uint Lo = Los[i];
		uint Hi = His[i];
		uint LoCol = MatchColToCol[Lo];
		uint HiCol = MatchColToCol[Hi];
		asserta(LoCol != UINT_MAX);
		asserta(HiCol != UINT_MAX);
		asserta(HiCol > LoCol);
		uint n = HiCol - LoCol + 1;
		Log("\n");
		Log("LSB[%u]\n", i);
		Log("%*.*s\n", n, n, RowQ.c_str() + LoCol);
		Log("%*.*s\n", n, n, RowT.c_str() + LoCol);
		}

	double Sum = 0;
	for (uint i = 0; i < LSBCount; ++i)
		{
		uint Loi = Los[i];
		uint Hii = His[i];
		for (uint j = 0; j < LSBCount; ++j)
			{
			uint Loj = Los[j];
			uint Hij = His[j];
			double Score = GetDALIScore_LSB_Pair(ProfQ, ProfT, PosQs, PosTs,
			  Loi, Hii, Loj, Hij);
			Sum += Score;
			Log("LSB[%u,%u] = %.4g  Sum %.4g\n", i, j, Score, Sum);
			}
		}
	Log("Sum = %.4g\n", Sum);

	uint Lo0 = Los[0];
	uint Hi0 = His[0];
	for (uint MatchCol = 0; MatchCol < MatchColCount; ++MatchCol)
		{
		double Score = GetDALIScore_LSB_Pair(ProfQ, ProfT, PosQs, PosTs,
		  Lo0, Hi0, MatchCol, MatchCol);
		Log("MatchCol %3u Score = %.4g\n", MatchCol, Score);
		}
	}
