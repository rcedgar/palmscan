#pragma once

#include "pdbchain.h"
#include "abcxyz.h"
#include <map>

const uint PPSPL = 34;
const uint AIX = 0;
const uint BIX = AL;
const uint CIX = AL + BL;

class CMP
	{
public:
	string m_Label;
	vector<string> m_RefLabels;
	vector<vector<double> > m_RefMeans;
	vector<vector<double> > m_StdDevs;
	vector<vector<vector<double> > > m_DistMxVec;

public:
	CMP()
		{
		Clear();
		};

	void Clear()
		{
		m_Label.clear();
		m_RefMeans.clear();
		m_StdDevs.clear();
		m_RefMeans.resize(PPSPL);
		m_StdDevs.resize(PPSPL);
		m_RefLabels.clear();
		for (uint i = 0; i < PPSPL; ++i)
			{
			m_RefMeans[i].resize(PPSPL, DBL_MAX);
			m_StdDevs[i].resize(PPSPL, DBL_MAX);
			}
		}

	void ToFile(const string &FileName) const;
	void MxToFile(FILE *f, const string &Name,
	  const vector<vector<double> > &Mx) const;
	void MxFromFile(FILE *f, string &Name,
	  vector<vector<double> > &Mx);
	void FromFile(FILE *f);
	void FromFile(const string &FileName);

// For training
public:
	void InitTrain()
		{
		m_DistMxVec.clear();
		}

	void AppendTrainDistMx(const vector<vector<double> > &DistMx)
		{
		asserta(SIZE(DistMx) == PPSPL);
		asserta(SIZE(DistMx[0]) == PPSPL);
		m_DistMxVec.push_back(DistMx);
		}

	void FinalizeTrain();
	void TrainChain(const PDBChain &Chain,
	  uint APos, uint BPos, uint CPos);
	void GetDistMx(const PDBChain &Chain, 
	  uint APos, uint BPos, uint CPos,
	  vector<vector<double> > &DistMx);
	void GetMeanStdDev(uint i, uint j,
	  double &Mean, double &StdDev) const;
	void ToDBFile(const string &FileName) const;

public:
	static char GetMotifChar(uint Ix);
	static uint GetSeqPos(uint i, uint APos, uint BPos, uint CPos);
	static double GetRdRpProb(int Gate, const string &GDD);
	static double GetRdRpProb_Gate(int Gate);
	static double GetRdRpProb_GDD(const string &GDD);
	};

double GetNormal(double Mu, double Sigma, double x);
void ReadMotifCoords(
  vector<vector<uint> > &MotifCoordsVec,
  vector<string> &Labels,
  map<string, uint> &LabelToIndex);
