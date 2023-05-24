#pragma once

#include "pdbchain.h"
#include "pf.h"

class PPP
	{
public:
	string m_Name;
	vector<vector<double> > m_FeatureValuesVec;
	vector<double> m_FeatureMeans;
	vector<double> m_FeatureStdDevs;

public:
	static uint m_ShapeLength;
	static uint m_ShapeSmoothWindow;
	static vector<string> m_FeatureNames;
	static uint m_FeatureCount;

public:
	void Clear()
		{
		m_Name.clear();
		m_FeatureValuesVec.clear();
		m_FeatureMeans.clear();
		m_FeatureStdDevs.clear();
		}

	void FromPPP(const string &FileName);
	double ScoreChain(const PDBChain &Chain,
	  uint PosA, uint PosB, uint PosC) const;
	void LogMe() const;
	void Init(const string &Name);
	void GetFeatureVec(const PDBChain &Chain,
	  vector<double> &FeatureValues) const;
	void AddChain(const PDBChain &Chain);
	void Averages();
	void Average(uint i);
	void ToPPP(FILE *f) const;
	void ToPDB(FILE *f) const;

public:
	static const char *GetFeatureName(uint i);
	static void GlobalInit();
	};
