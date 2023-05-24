#include "myutils.h"
#include "outputfiles.h"
#include "sort.h"
#include "quarts.h"
#include "xtrainer.h"
#include "pdbchain.h"
#include "xprof.h"

void XTrainer::LogProbs() const
	{
	Log("\n");
	Log("____________ %s ____________\n",
	  XProf::GetFeatureName(m_FeatureIndex));
	for (uint i = 0; i < XBINS; ++i)
		Log("[%u]  %6.1f%%\n", i, 100.0*m_Probs[i]);

	Log("\n");
	for (uint i = 0; i < XBINS; ++i)
		{
		Log("[%u]", i);
		for (uint j = 0; j < XBINS; ++j)
			{
			double P = m_JointProbs[i][j];
			asserta(m_JointProbs[i][j] == P);
			Log("  %6.2g", P);
			}
		Log("\n");
		}

	Log("\n");
	for (uint i = 0; i < XBINS; ++i)
		{
		Log("[%u]", i);
		for (uint j = 0; j < XBINS; ++j)
			{
			double Score = m_LogOddsMx[i][j];
			asserta(m_LogOddsMx[j][i] == Score);
			Log("  %4.1f", Score);
			}
		Log("\n");
		}
	}

void XTrainer::LogFinal() const
	{
	const uint FeatureCount = XProf::GetFeatureCount();
	
	Log("\n");
	Log("Feature bin pcts.\n");
	for (uint Bin = 0; Bin < XBINS; ++Bin)
		Log("%7u", Bin);
	Log("\n");
	for (uint Bin = 0; Bin < XBINS; ++Bin)
		Log("%7.7s", "-----");
	Log("\n");
	for (uint i = 0; i < FeatureCount; ++i)
		{
		const char *Name = XProf::GetFeatureName(i);
		for (uint Bin = 0; Bin < XBINS; ++Bin)
			{
			double Prob = m_FeatureIndexToProbs[i][Bin];
			Log("%7.1f", 100.0*Prob);
			}
		Log("  |  %s\n", Name);
		}

	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		const vector<vector<double> > &LogOddsMx = m_FeatureIndexToLogOddsMx[FeatureIndex];
		Log("\n");
		Log("_________ %s _________\n", XProf::GetFeatureName(FeatureIndex));
		for (uint i = 0; i < XBINS; ++i)
			{
			Log("[%u]", i);
			for (uint j = 0; j <= i; ++j)
				{
				double Score = LogOddsMx[i][j];
				asserta(LogOddsMx[j][i] == Score);
				Log("  %4.1f", Score);
				}
			Log("\n");
			}
		}
	}

uint XTrainer::GetSeqLength(uint ChainIndex) const
	{
	asserta(ChainIndex < SIZE(m_Chains));
	return m_Chains[ChainIndex]->GetSeqLength();
	}

const string &XTrainer::GetLabel(uint ChainIndex) const
	{
	asserta(ChainIndex < SIZE(m_Chains));
	return m_Chains[ChainIndex]->m_Label;
	}

uint XTrainer::GetIndex(const string &Label) const
	{
	map<string, uint>::const_iterator p =
	  m_LabelToChainIndex.find(Label);
	if (p == m_LabelToChainIndex.end())
		return UINT_MAX;
	return p->second;
	}

void XTrainer::SetChains(const vector<PDBChain *> &Chains)
	{
	m_Chains.clear();
	m_XProfs.clear();

	vector<string> Labels;
	uint N = GetChainCount();
	for (uint i = 0; i < SIZE(Chains); ++i)
		{
		PDBChain *Chain = Chains[i];
		XProf *XP = new XProf;
		XP->Init(*Chain);
		m_Chains.push_back(Chain);
		m_XProfs.push_back(XP);
		const string &Label = Chain->m_Label;
		m_LabelToChainIndex[Label] = i;
		Labels.push_back(Label);
		}
	}

void XTrainer::SetAln(const string &Label1, const string &Label2,
  const string &Row1, const string &Row2)
	{
	uint ChainIndex1 = GetIndex(Label1);
	uint ChainIndex2 = GetIndex(Label2);
	if (ChainIndex1 == UINT_MAX || ChainIndex2 == UINT_MAX)
		return;

	string S1 = Row1;
	string S2 = Row2;

	StripGaps(S1);
	StripGaps(S2);

	if (S1 != m_Chains[ChainIndex1]->m_Seq)
		{
		Log("\n");
		Log("Aln %s - %s\n", Label1.c_str(), Label2.c_str());
		Log("%s\n", S1.c_str());
		Log("%s\n", m_Chains[ChainIndex1]->m_Seq.c_str());
		Warning("S1");
		}

	if (S2 != m_Chains[ChainIndex2]->m_Seq)
		{
		Log("\n");
		Log("Aln %s - %s\n", Label1.c_str(), Label2.c_str());
		Log("%s\n", S2.c_str());
		Log("%s\n", m_Chains[ChainIndex2]->m_Seq.c_str());
		Warning("S2");
		}

	const pair<uint, uint> Pair12(ChainIndex1, ChainIndex2);
	const pair<string, string> RowPair12(Row1, Row2);
	m_PairToAlnRows[Pair12] = RowPair12;

	const pair<uint, uint> Pair21(ChainIndex2, ChainIndex1);
	const pair<string, string> RowPair21(Row2, Row1);
	m_PairToAlnRows[Pair21] = RowPair21;
	}

bool XTrainer::GetAln(uint ChainIndex1, uint ChainIndex2,
  string &Row1, string &Row2) const
	{
	Row1.clear();
	Row2.clear();
	const pair<uint, uint> Pair(ChainIndex1, ChainIndex2);
	map<pair<uint, uint>, pair<string, string> >::const_iterator p =
	  m_PairToAlnRows.find(Pair);
	if (p == m_PairToAlnRows.end())
		return false;
	Row1 = p->second.first;
	Row2 = p->second.second;
	return true;
	}

void XTrainer::InitTrainFeature(uint FeatureIndex)
	{
	m_FeatureIndex = FeatureIndex;

	m_Counts.clear();
	m_JointCounts.clear();
	m_Probs.clear();
	m_JointProbs.clear();
	m_LogOddsMx.clear();

	m_Probs.resize(XBINS, 0);
	m_Counts.resize(XBINS, 0);

	m_JointCounts.resize(XBINS);
	m_JointProbs.resize(XBINS);
	m_LogOddsMx.resize(XBINS);
	for (uint i = 0; i < XBINS; ++i)
		{
		m_JointCounts[i].resize(XBINS, 0);
		m_LogOddsMx[i].resize(XBINS, DBL_MAX);
		m_JointProbs[i].resize(XBINS, DBL_MAX);
		}
	}

void XTrainer::CalcProbsFeature()
	{
	uint Sum = 0;
	for (uint i = 0; i < XBINS; ++i)
		Sum += m_Counts[i];
	asserta(Sum > 0);

	uint Sum2 = 0;
	for (uint i = 0; i < XBINS; ++i)
		for (uint j = 0; j < XBINS; ++j)
			Sum2 += m_JointCounts[i][j];
	asserta(Sum == Sum2);

	for (uint i = 0; i < XBINS; ++i)
		{
		uint n = m_Counts[i];
		double P = double((n+1))/(Sum+1);
		m_Probs[i] = P;
		}

	for (uint i = 0; i < XBINS; ++i)
		{
		double Pi = m_Probs[i];
		for (uint j = 0; j < XBINS; ++j)
			{
			double Pj = m_Probs[j];
			double n = (m_JointCounts[i][j] + m_JointCounts[j][i])/2.0;
			double P = double(n)/Sum;
			double r = P/(Pi*Pj);
			double Score = 0;
			if (r < 0.001)
				Score = -3;
			else if (r > 100)
				Score = +3;
			else
				Score = 10*log10(r);

			m_JointProbs[i][j] = P;
			//m_JointProbs[j][i] = P;
			m_LogOddsMx[i][j] = Score;
			//m_LogOddsMx[j][i] = Score;
			}
		}
	}

void XTrainer::TrainFeature(uint FeatureIndex)
	{
	InitTrainFeature(FeatureIndex);

	const uint N = GetChainCount();
	for (uint i = 1; i < N; ++i)
		for (uint j = 0; j < i; ++j)
			TrainFeaturePair(i, j);

	CalcProbsFeature();
	}

void XTrainer::TrainFeaturePair(uint ChainIndex1, uint ChainIndex2)
	{
	string Row1, Row2;
	bool Ok = GetAln(ChainIndex1, ChainIndex2, Row1, Row2);
	if (!Ok)
		return;
	uint ColCount = SIZE(Row1);
	asserta(SIZE(Row2) == ColCount);

	const uint L1 = GetSeqLength(ChainIndex1);
	const uint L2 = GetSeqLength(ChainIndex2);

	uint Pos1 = 0;
	uint Pos2 = 0;

	for (uint Col = 0; Col < ColCount; ++Col)
		{
		char c1 = Row1[Col];
		char c2 = Row2[Col];

		if (c1 != '-' && c2 != '-')
			{
			TrainFeaturePairPos(ChainIndex1, ChainIndex2,
			  Pos1, Pos2, true);
			TrainFeaturePairPos(ChainIndex1, ChainIndex2,
			  Pos1, L2 - Pos2 - 1, false);
			}

		if (c1 != '-')
			++Pos1;
		if (c2 != '-')
			++Pos2;
		}
	asserta(Pos1 == L1);
	asserta(Pos2 == L2);
	}

void XTrainer::TrainFeaturePairPos(uint ChainIndex1, uint ChainIndex2,
  uint Pos1, uint Pos2, bool TOF)
	{
	const XProf &XP1 = *m_XProfs[ChainIndex1];
	const XProf &XP2 = *m_XProfs[ChainIndex2];
	const uint L1 = GetSeqLength(ChainIndex1);
	const uint L2 = GetSeqLength(ChainIndex2);
	asserta(Pos1 < L1);
	asserta(Pos2 < L2);
	double v1;
	double v2;
	uint i1;
	uint i2;
	XP1.GetFeature(m_FeatureIndex, Pos1, v1, i1);
	XP2.GetFeature(m_FeatureIndex, Pos2, v2, i2);
	if (i1 == UINT_MAX || i2 == UINT_MAX)
		return;

	asserta(i1 < XBINS);
	asserta(i2 < XBINS);
	++(m_Counts[i1]);
	++(m_JointCounts[i1][i2]);
	}

void XTrainer::SetValues(uint FeatureIndex)
	{
	m_FeatureIndex = FeatureIndex;
	m_Values.clear();
	m_iValues.clear();
	const uint N = GetChainCount();
	for (uint ChainIndex = 0; ChainIndex < N; ++ChainIndex)
		SetValues1(ChainIndex);
	}

void XTrainer::LogValueHist() const
	{
	Quarts Q;
	QuartsDouble QD;
	GetQuarts(m_iValues, Q);
	GetQuartsDouble(m_Values, QD);
	Q.LogMe();
	QD.LogMe();
	}

void XTrainer::SetValues1(uint ChainIndex)
	{
	const uint L = GetSeqLength(ChainIndex);
	const XProf &XP = *m_XProfs[ChainIndex];
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		double Value;
		uint iValue;
		XP.GetFeature(m_FeatureIndex, Pos, Value, iValue);
		if (Value != DBL_MAX)
			{
			asserta(Value >= 0 && Value < 999);
			m_Values.push_back(Value);
			}
		if (iValue != UINT_MAX)
			{
			asserta(iValue < XBINS);
			m_iValues.push_back(iValue);
			}
		}
	}

void XTrainer::Train()
	{
	m_FeatureIndexToValues.clear();
	m_FeatureIndexToProbs.clear();
	m_FeatureIndexToJointProbs.clear();
	m_FeatureIndexToLogOddsMx.clear();

	const uint FeatureCount = XProf::GetFeatureCount();
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
	  ++FeatureIndex)
		{
		ProgressStep(FeatureIndex, FeatureCount, "Train");
		TrainFeature(FeatureIndex);

//		LogProbs();

		m_FeatureIndexToValues.push_back(m_Values);
		m_FeatureIndexToProbs.push_back(m_Probs);
		m_FeatureIndexToJointProbs.push_back(m_JointProbs);
		m_FeatureIndexToLogOddsMx.push_back(m_LogOddsMx);
		}
	}

double XTrainer::CorrelFeatures(uint FeatureIndex1, uint FeatureIndex2) const
	{
	double SumX = 0;
	double SumY = 0;
	double SumX2 = 0;
	double SumY2 = 0;
	double SumXY = 0;

	const uint ChainCount = GetChainCount();
	uint N = 0;
	for (uint i = 0; i < ChainCount; ++i)
		{
		const XProf &XP = *m_XProfs[i];
		const uint L = GetSeqLength(i);
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			double Value1;
			double Value2;
			uint iValue1;
			uint iValue2;
			XP.GetFeature(FeatureIndex1, Pos, Value1, iValue1);
			XP.GetFeature(FeatureIndex2, Pos, Value2, iValue2);
			if (Value1 == DBL_MAX || Value2 == DBL_MAX)
				continue;

			++N;
			SumX += Value1;
			SumY += Value2;
			SumX2 += Value1*Value1;
			SumY2 += Value2*Value2;
			SumXY += Value1*Value2;
			}
		}
	asserta(N > 0);
	double Top = SumXY - (SumX*SumY)/N;
	double BottomL = SumX2 - (SumX*SumX)/N;
	double BottomR = SumY2 - (SumY*SumY)/N;
	double r = Top/sqrt(BottomL*BottomR);
	return r;
	}

double XTrainer::GetScorePosPair(uint ChainIndex1, uint Pos1, 
  uint ChainIndex2, uint Pos2) const
	{
	const uint FeatureCount = XProf::GetFeatureCount();
	const XProf &XP1 = *m_XProfs[ChainIndex1];
	const XProf &XP2 = *m_XProfs[ChainIndex2];

	double Sum = 0;
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount; ++FeatureIndex)
		{
		double Value1;
		double Value2;
		uint iValue1;
		uint iValue2;
		XP1.GetFeature(FeatureIndex, Pos1, Value1, iValue1);
		XP2.GetFeature(FeatureIndex, Pos2, Value2, iValue2);
		if (iValue1 != UINT_MAX && iValue2 != UINT_MAX)
			{
			double Score =
			  m_FeatureIndexToLogOddsMx[FeatureIndex][iValue1][iValue2];
			Sum += Score;
			}
		}
	return Sum/FeatureCount;
	}

void XTrainer::GetScoreMx(uint ChainIndex1, uint ChainIndex2,
  vector<vector<double> > &Mx) const
	{
	Mx.clear();
	const PDBChain &Chain1 = *m_Chains[ChainIndex1];
	const PDBChain &Chain2 = *m_Chains[ChainIndex2];
	const uint L1 = Chain1.GetSeqLength();
	const uint L2 = Chain2.GetSeqLength();

	Mx.resize(L1);
	for (uint Pos1 = 0; Pos1 < L1; ++Pos1)
		{
		for (uint Pos2 = 0; Pos2 < L2; ++Pos2)
			{
			double Score =
			  GetScorePosPair(ChainIndex1, Pos1, ChainIndex2, Pos2);
			Mx[Pos1].push_back(Score);
			}
		asserta(SIZE(Mx[Pos1]) == L2);
		}
	}

void XTrainer::GetScoreMx2(uint ChainIndex1, uint ChainIndex2,
  Mx<float> &SMx) const
	{
	const PDBChain &Chain1 = *m_Chains[ChainIndex1];
	const PDBChain &Chain2 = *m_Chains[ChainIndex2];
	const uint L1 = Chain1.GetSeqLength();
	const uint L2 = Chain2.GetSeqLength();
	SMx.Alloc("SMx", L1, L2);

	for (uint Pos1 = 0; Pos1 < L1; ++Pos1)
		{
		for (uint Pos2 = 0; Pos2 < L2; ++Pos2)
			{
			float Score = (float)
			  GetScorePosPair(ChainIndex1, Pos1, ChainIndex2, Pos2);
			SMx.Put(Pos1, Pos2, Score);
			}
		}
	}

void XTrainer::LogOddsToTsv(const string &FileName,
  vector<vector<vector<double> > > &FeatureIndexToLogOddsMx)
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	const uint FeatureCount = XProf::GetFeatureCount();
	fprintf(f, "logodds\t%u\n", FeatureCount);
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
	   ++FeatureIndex)
		{
		const char *Name = XProf::GetFeatureName(FeatureIndex);
	
		fprintf(f, "logodds1\t%u\t%u\t%s\n",
		  FeatureIndex, XBINS, Name);

		asserta(FeatureIndex < SIZE(FeatureIndexToLogOddsMx));
		const vector<vector<double> > &LogOddsMx =
		  FeatureIndexToLogOddsMx[FeatureIndex];
		asserta(SIZE(LogOddsMx) == XBINS);
		for (uint i = 0; i < XBINS; ++i)
			{
			fprintf(f, "%u", i);
			const vector<double> &Row = LogOddsMx[i];
			for (uint j = 0; j <= i; ++j)
				fprintf(f, "\t%.4g", Row[j]);
			fprintf(f, "\n");
			}
		}
	}

void XTrainer::LogOddsFromTsv(const string &FileName,
  vector<vector<vector<double> > > &FeatureIndexToLogOddsMx)
	{
	FeatureIndexToLogOddsMx.clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	asserta(Fields[0] == "logodds");
	uint FeatureCount = StrToUint(Fields[1]);
	asserta(FeatureCount == XProf::GetFeatureCount());

	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
	   ++FeatureIndex)
		{
		const string &Name = XProf::GetFeatureName(FeatureIndex);

		Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 4);
		asserta(Fields[0] == "logodds1");
		asserta(StrToUint(Fields[1]) == FeatureIndex);
		uint xb = StrToUint(Fields[2]);
		asserta(xb == XBINS);
		asserta(Fields[3] == Name);

		vector<vector<double> > LogOddsMx(xb);
		for (uint i = 0; i < xb; ++i)
			LogOddsMx[i].resize(xb, DBL_MAX);

		for (uint i = 0; i < xb; ++i)
			{
			Ok = ReadLineStdioFile(f, Line);
			asserta(Ok);
			Split(Line, Fields, '\t');
			asserta(SIZE(Fields) == i+2);
			asserta(StrToUint(Fields[0]) == i);
			for (uint j = 0; j <= i; ++j)
				{
				double Score = StrToFloat(Fields[j+1]);
				LogOddsMx[i][j] = Score;
				LogOddsMx[j][i] = Score;
				}
			}
		FeatureIndexToLogOddsMx.push_back(LogOddsMx);
		}
	}
