#pragma once

#include "pdbchain.h"

class StructProf
	{
public:
	const PDBChain *m_Chain = 0;
	//uint m_MinPos = UINT_MAX;
	//uint m_MaxPos = UINT_MAX;
	vector<double> m_CavityCenterPt;

public:
	void Clear()
		{
		m_Chain = 0;
		//m_MinPos = UINT_MAX;
		//m_MaxPos = UINT_MAX;
		m_CavityCenterPt.clear();
		}

	void SetChain(const PDBChain &Chain);
	//void SetMinMaxPos(uint MinPos, uint MaxPos);
	void SetCavityCenterPt();
	void GetHSE(uint Pos, double Radius,
	  uint &NU, uint &ND) const;
	uint GetTSB(uint Pos, double Radius) const;
	uint GetCavityNumber(uint Pos) const;
	void GetSphere(const vector<double> &CenterPt, double Radius,
	  vector<uint> &PosVec) const;
	uint SearchDist(uint Pos, uint Lo, uint Hi,
	  bool Maximize, double X, double &BestDist,
	  bool Trace = false) const;
	uint FindMofifD_Hueuristics() const;
	uint FindMofifE_Hueuristics(uint Pos_MotifD) const;
	uint FindMofifF1_Hueuristics(uint Pos_MotifA) const;
	uint FindMofifF2_Hueuristics(uint Pos_MotifA) const;
	void WriteGSProf(FILE *f) const;
	void WriteTsv(FILE *f) const;
	void WriteMotifTsv(FILE *f, uint Pos, uint n) const;
	};
