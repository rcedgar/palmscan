#pragma once

#include <map>
#include "mx.h"

class XProf;
class PDBChain;

class XTrainer
	{
public:
	vector<PDBChain *> m_Chains;
	vector<XProf *> m_XProfs;
	map<string, uint> m_LabelToChainIndex;
	map<pair<uint, uint>, pair<string, string> > m_PairToAlnRows;

	uint m_FeatureIndex = UINT_MAX;
	vector<double> m_Values;
	vector<uint> m_iValues;
	vector<uint> m_Counts;
	vector<vector<uint> > m_JointCounts;
	vector<double> m_Probs;
	vector<vector<double> > m_JointProbs;
	vector<vector<double> > m_LogOddsMx;

	vector<vector<double> > m_FeatureIndexToValues;
	vector<vector<double> > m_FeatureIndexToProbs;
	vector<vector<vector<double> > > m_FeatureIndexToJointProbs;
	vector<vector<vector<double> > > m_FeatureIndexToLogOddsMx;

public:
	void LogValueHist() const;
	void LogProbs() const;
	void LogFinal() const;
	void InitTrainFeature(uint FeatureIndex);
	void CalcProbsFeature();
	uint GetChainCount() const { return SIZE(m_Chains); }
	void SetChains(const vector<PDBChain *> &Chains);
	uint GetIndex(const string &Label) const;
	void SetAln(const string &Label1, const string &Label2,
	  const string &Row1, const string &Row2);
	bool GetAln(uint ChainIndex1, uint ChainIndex2, string &Row1, string &Row2) const;
	uint GetSeqLength(uint ChainIndex) const;
	const string &GetLabel(uint ChainIndex) const;
	void SetValues(uint FeatureIndex);
	void SetValues1(uint ChainIndex);
	void Train();
	void TrainFeature(uint FeatureIndex);
	void TrainFeaturePair(uint ChainIndex1, uint ChainIndex2);
	void TrainFeaturePairPos(uint ChainIndex1, uint ChainIndex2, uint Pos1, uint Pos2, bool TOF);
	double CorrelFeatures(uint FeatureIndex1, uint FeatureIndex2) const;
	void GetScoreMx(uint ChainIndex1, uint ChainIndex2,
	  vector<vector<double> > &Mx) const;
	void GetScoreMx2(uint ChainIndex1, uint ChainIndex2,
	  Mx<float> &SMx) const;
	double GetScorePosPair(uint ChainIndex1, uint Pos1, 
	  uint ChainIndex2, uint Pos2) const;

public:
	static void LogOddsToTsv(const string &FileName,
	  vector<vector<vector<double> > > &FeatureIndexToLogOddsMx);
	static void LogOddsFromTsv(const string &FileName,
	  vector<vector<vector<double> > > &FeatureIndexToLogOddsMx);
	};
