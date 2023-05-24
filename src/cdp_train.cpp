#include "myutils.h"
#include "abcxyz.h"
#include "outputfiles.h"
#include "cddata.h"
#include "cdsearcher.h"
#include "quarts.h"
#include <map>

void SetPalmTemplate(const CDInfo &Info, CDTemplate &Tpl)
	{
	Tpl.Clear();

	vector<char> AnchorAAs;
	vector<uint> AnchorAAOffsets;
	vector<uint> MinAAsNext;
	vector<uint> MaxAAsNext;
	vector<double> MinScores;

/***
On palmcore training set:
F2  minaa  51  maxaa  74  minscore 0.8094
 A  minaa  52  maxaa 130  minscore 0.6463
 B  minaa  30  maxaa  75  minscore 0.7778
 C  minaa  24  maxaa  45  minscore 0.7640
 D  minaa   8  maxaa  21  minscore 0.7365
 E  minaa   -  maxaa   -  minscore 0.5780
***/
#define S(AA, Off, MinNext, MaxNext, MinScore) \
	AnchorAAs.push_back(AA); \
	AnchorAAOffsets.push_back(Off); \
	MinAAsNext.push_back(MinNext); \
	MaxAAsNext.push_back(MaxNext); \
	MinScores.push_back(MinScore);

	S('R', 2, 30, 100, 0.6);
	S('D', 3, 30, 150, 0.55);
	S('G', 1, 20, 100, 0.6);
	S('D', 3, 10, 70, 0.6);
	S('x', 0, 4, 70, 0.6);
	S('x', 0, 0, 100, 0.4);

#undef S

	Tpl.Init(Info, AnchorAAs, AnchorAAOffsets,
	  MinAAsNext, MaxAAsNext, MinScores);
	}

static void ReadMotifDataCoords(
  CDInfo &Info,
  vector<vector<uint> > &MotifCoordsVec,
  vector<vector<string> > &MotifSeqsVec,
  vector<string> &Labels,
  map<string, uint> &LabelToIndex)
	{
	MotifCoordsVec.clear();
	MotifSeqsVec.clear();
	Labels.clear();
	LabelToIndex.clear();

	vector<string> MotifNames;
	vector<uint> MotifLengths;

	MotifNames.push_back("F2");
	MotifNames.push_back("A");
	MotifNames.push_back("B");
	MotifNames.push_back("C");
	MotifNames.push_back("D");
	MotifNames.push_back("E");

	MotifLengths.push_back(7);
	MotifLengths.push_back(12);
	MotifLengths.push_back(14);
	MotifLengths.push_back(10);
	MotifLengths.push_back(7);
	MotifLengths.push_back(7);

	const uint NM = 6;

	Info.Init(MotifNames, MotifLengths);

	if (!optset_motif_coords)
		Die("Must specify -motif_coords");

	FILE *f = OpenStdioFile(opt_motif_coords);
	string Line;
	vector<string> Fields;
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 2*NM + 1);

		vector<uint> MotifCoords;
		vector<string> MotifSeqs;

		const string &Label = Fields[0];

		for (uint MotifIndex = 0; MotifIndex < NM; ++MotifIndex)
			{
			const string &s = Fields[1 + 2*MotifIndex];
			if (s != ".")
				{
				uint Pos = StrToUint(s);
				asserta(Pos > 0);
				MotifCoords.push_back(Pos-1);
				}
			else
				MotifCoords.push_back(UINT_MAX);
			MotifSeqs.push_back(Fields[2 + 2*MotifIndex]);
			}

		uint Index = SIZE(Labels);
		asserta(LabelToIndex.find(Label) == LabelToIndex.end());
		LabelToIndex[Label] = Index;
		Labels.push_back(Label);
		MotifCoordsVec.push_back(MotifCoords);
		MotifSeqsVec.push_back(MotifSeqs);
		}
	}

static void GetData(const PDBChain &Chain,
  const CDInfo &Info,
  const vector<uint> &MotifCoords, 
  const vector<string> &MotifSeqs,
  CDData &Data)
	{
	//Log("\n_________________________\n");
	//Log(">%s\n", Chain.m_Label.c_str());

	const uint NM = Info.m_MotifCount;
	for (uint MotifIndex = 0; MotifIndex < NM; ++MotifIndex)
		{
		uint ML = Info.m_MotifLengths[MotifIndex];
		uint Pos = MotifCoords[MotifIndex];
		if (Pos == UINT_MAX)
			continue;
		const string &MotifSeq = MotifSeqs[MotifIndex];
		string Motif;
		Chain.GetSubSeq(Pos, ML, Motif);
		//Log(" %3.3s", Info.m_MotifNames[MotifIndex].c_str());
		//Log(" %4u", Pos + 1);
		//Log("  %16.16s", MotifSeq.c_str());
		//Log("  %16.16s", Motif.c_str());
		//Log("\n");
		asserta(Motif == MotifSeq);
		}

	Data.Init(Info);
	const uint Size = Info.GetSize();
	for (uint Ix1 = 0; Ix1 < Size; ++Ix1)
		{
		uint MotifIndex1 = Info.m_MotifIndexes[Ix1];
		asserta(MotifIndex1 < NM);
		uint Coord1 = Info.m_MotifCoords[Ix1];
		uint MotifPos1 = MotifCoords[MotifIndex1];
		if (MotifPos1 == UINT_MAX)
			continue;
		uint Pos1 = MotifPos1 + Coord1;
		for (uint Ix2 = 0; Ix2 < Size; ++Ix2)
			{
			uint MotifIndex2 = Info.m_MotifIndexes[Ix2];
			asserta(MotifIndex2 < NM);
			uint Coord2 = Info.m_MotifCoords[Ix2];
			uint MotifPos2 = MotifCoords[MotifIndex2];
			if (MotifPos2 == UINT_MAX)
				continue;
			uint Pos2 = MotifPos2 + Coord2;
			double Dist = Chain.GetDist(Pos1, Pos2);
			Data.SetByIx(Ix1, Ix2, Dist);
			}
		}
	}

static void GetAvgStdDev(const vector<CDData> &DataVec,
  const CDInfo &Info, CDData &DataAvg, CDData &DataStdDev)
	{
	DataAvg.Init(Info);
	DataStdDev.Init(Info);

	const uint N = SIZE(DataVec);
	asserta(N > 0);
	const uint Size = Info.GetSize();
	for (uint Ix1 = 0; Ix1 < Size; ++Ix1)
		{
		for (uint Ix2 = 0; Ix2 < Size; ++Ix2)
			{
			vector<double> Dists;
			for (uint i = 0; i < N; ++i)
				{
				double Dist = DataVec[i].GetByIx(Ix1, Ix2);
				if (Dist == DBL_MAX)
					continue;
				Dists.push_back(Dist);
				}

			asserta(!Dists.empty());
			QuartsDouble QD;
			GetQuartsDouble(Dists, QD);
			double Avg = QD.Avg;
			double StdDev = QD.StdDev;
			DataAvg.SetByIx(Ix1, Ix2, Avg);
			DataStdDev.SetByIx(Ix1, Ix2, StdDev);
			}
		}
	}

static uint GetBestFit(const CDSearcher &CS, uint MotifIndex,
  uint TrainPos, double &BestScore)
	{
	BestScore = 0;
	if (TrainPos == UINT_MAX)
		return UINT_MAX;
	const int iQL = (int) CS.m_Query->GetSeqLength();
	const CDInfo &Info = *CS.m_Info;
	int ML = Info.GetMotifLength(MotifIndex);
	uint Ix = Info.GetIx(MotifIndex, 0);
	int iTrainPos = int(TrainPos);
	int iBestPos = -999;
	for (int iPos2 = iTrainPos - 4; iPos2 < iTrainPos + 4; ++iPos2)
		{
		if (iPos2 < 0 || iPos2 + (int) ML >= iQL)
			continue;
		uint Pos2 = uint(iPos2);
		double Score = CS.GetScore(Pos2, Pos2, Ix, Ix, ML, ML);
		if (Score > BestScore)
			{
			iBestPos = iPos2;
			BestScore = Score;
			}
		}
	if (iBestPos == -999)
		return UINT_MAX;
	return (uint) iBestPos;
	}

static void UpdateMinMax(const vector<uint> &MotifCoords, 
  vector<uint> &MotifIndexToMinAAs,
  vector<uint> &MotifIndexToMaxAAs)
	{
	uint MotifCount = SIZE(MotifCoords);
	asserta(SIZE(MotifIndexToMinAAs) == MotifCount);
	asserta(SIZE(MotifIndexToMaxAAs) == MotifCount);
	for (uint MotifIndex = 0; MotifIndex + 1 < MotifCount;
	  ++MotifIndex)
		{
		uint Pos = MotifCoords[MotifIndex];
		uint PosNext = MotifCoords[MotifIndex+1];
		if (Pos == UINT_MAX || PosNext == UINT_MAX)
			continue;

		asserta(PosNext > Pos);
		uint AAs = PosNext - Pos;

		MotifIndexToMinAAs[MotifIndex] =
		  min(AAs, MotifIndexToMinAAs[MotifIndex]);

		MotifIndexToMaxAAs[MotifIndex] =
		  max(AAs, MotifIndexToMaxAAs[MotifIndex]);
		}
	}

void LogCDHit(const string &Msg, const CDSearcher &CS,
  const PDBChain &Q, const vector<uint> &MotifCoords)
	{
	Log("%16.16s", Q.m_Label.c_str());
	Log("  %8.8s", Msg.c_str());
	if (MotifCoords.empty())
		{
		Log("  (empty)\n");
		return;
		}
	const uint MotifCount = CS.m_Info->GetMotifCount();
	asserta(SIZE(MotifCoords) == MotifCount);
	vector<uint> MotifIndexes;
	bool All = true;
	for (uint MotifIndex = 0; MotifIndex < MotifCount; ++MotifIndex)
		{
		MotifIndexes.push_back(MotifIndex);
		uint ML = CS.m_Info->GetMotifLength(MotifIndex);
		uint Pos = MotifCoords[MotifIndex];
		if (Pos == UINT_MAX)
			{
			All = false;
			Log("  %4.4s  %*.*s", "", ML, ML, "");
			continue;
			}

		string Seq;
		Q.GetSubSeq(Pos, ML, Seq);
		Log("  %4u  %*.*s", Pos+1, ML, ML, Seq.c_str());
		}
	if (All)
		{
		double Score = CS.GetScoreHit(MotifIndexes, MotifCoords);
		Log("  %6.4f", Score);
		}
	Log("\n");
	}

void cmd_cdp_train()
	{
	const string &InputFN = opt_cdp_train;

	CDInfo Info;
	vector<vector<string> > MotifSeqsVec;
	vector<vector<uint> > MotifCoordsVec;
	vector<string> Labels;
	map<string, uint> LabelToIndex;
	ReadMotifDataCoords(Info, MotifCoordsVec, 
	  MotifSeqsVec, Labels, LabelToIndex);
	const uint MotifCount = Info.GetMotifCount();

	PDBChain Q;
	CDData Data;
	vector<CDData> DataVec;

	vector<PDBChain *> Chains;
	ReadChains(InputFN, Chains);
	const uint N = SIZE(Chains);

	vector<uint> MotifIndexToMinAAs(MotifCount, UINT_MAX);
	vector<uint> MotifIndexToMaxAAs(MotifCount, 0);

//////////////////////////////////////////////////////
//  Train
//////////////////////////////////////////////////////
	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Q = *Chains[i];
		const string &Label = Q.m_Label;
		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			{
			Log("Not found >%s\n", Label.c_str());
			continue;
			}

		uint Index = p->second;
		const vector<uint> &MotifCoords = MotifCoordsVec[Index];
		const vector<string> &MotifSeqs = MotifSeqsVec[Index];

		UpdateMinMax(MotifCoords, MotifIndexToMinAAs, MotifIndexToMaxAAs);

		GetData(Q, Info, MotifCoords, MotifSeqs, Data);
		DataVec.push_back(Data);
		}

	CDData DataAvg;
	CDData DataStdDev;
	GetAvgStdDev(DataVec, Info, DataAvg, DataStdDev);

	if (g_ftsv != 0)
		{
		Info.ToTsv(g_ftsv);
		DataAvg.ToTsv("avg", g_ftsv);
		DataStdDev.ToTsv("stddev", g_ftsv);
		}

	vector<uint> MotifsABC;
	MotifsABC.push_back(1);
	MotifsABC.push_back(2);
	MotifsABC.push_back(3);

//////////////////////////////////////////////////////
//  Test
//////////////////////////////////////////////////////
	CDSearcher CS;
	CS.Init(Info, DataAvg, DataStdDev);

	CDTemplate Tpl;
	SetPalmTemplate(Info, Tpl);
	CS.m_Template = &Tpl;

	for (uint i = 0; i < N; ++i)
		{
		const PDBChain &Q = *Chains[i];
		CS.m_Query = &Q;
		const string &Label = Q.m_Label;
		const uint QL = Q.GetSeqLength();
		map<string, uint>::const_iterator p = LabelToIndex.find(Label);
		if (p == LabelToIndex.end())
			continue;

		uint Index = p->second;
		const vector<uint> &MotifCoords = MotifCoordsVec[Index];
		const vector<string> &MotifSeqs = MotifSeqsVec[Index];

		vector<uint> PalmHit;
		CS.SearchPalm(Q, PalmHit);

		Log("\n");
		LogCDHit("Train", CS, Q, MotifCoords);
		LogCDHit("Test", CS, Q, PalmHit);
		const char *Result = "Agree";
		if (PalmHit.empty())
			Result = "(train negative)";
		else
			{
			for (uint i = 0; i < MotifCount; ++i)
				{
				uint p1 = MotifCoords[i];
				uint p2 = PalmHit[i];
				if (p1 == UINT_MAX || p2 == UINT_MAX)
					{
					Result = "(missing motif)";
					break;
					}
				else if (p1 != p2)
					{
					Result = "DIFF";
					break;
					}
				}
			}
		Log(" -- %s\n", Result);
		}
	}
