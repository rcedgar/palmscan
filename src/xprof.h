#ifndef xprof_h
#define xprof_h

#include "myutils.h"
#include "pdbchain.h"

const uint XBINS = 8;

class XProf
	{
public:
	const PDBChain *m_Chain = 0;
	uint m_L = 0;

public:
	void Init(const PDBChain &Chain);
	void ToCfv(FILE *f) const;
	void PosToCfv(FILE *f, uint Pos) const;

	uint GetSeqLength() const { return m_Chain->GetSeqLength(); }

	void GetFeature(uint FeatureIndex, uint Pos,
	  double &Value, uint &iValue) const;

	void Get_NU(uint Pos, double &Value, uint &iValue) const;
	void Get_ND(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang_m1_p1(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang_m2_p2(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang_m3_p3(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang01_23(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang01_34(uint Pos, double &Value, uint &iValue) const;
	void Get_Ang01_45(uint Pos, double &Value, uint &iValue) const;
	void Get_ED_p4(uint Pos, double &Value, uint &iValue) const;
	void Get_ED_m4(uint Pos, double &Value, uint &iValue) const;

	uint IntScale(double Value, double MinVal, double MedVal) const;

	double GetAngle(
	  uint PosA1, uint PosA2, 
	  uint PosB1, uint PosB2) const;
	uint GetSphereNr(uint Pos, double Radius) const;

public:
	static uint GetFeatureCount();
	static const char *GetFeatureName(uint FeatureIndex);
	};

#endif // xprof_h
