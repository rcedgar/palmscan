#include "myutils.h"
#include "chainreader.h"
#include "outputfiles.h"
#include "abcxyz.h"
#include "pdbchain.h"
#include "xprof.h"

void GetHSE(const PDBChain &Chain, uint Pos, double Radius,
  uint MinPos, uint MaxPos, uint &NU, uint &ND);

void cmd_xprof()
	{
	const string &InputFN = opt_xprof;

	ChainReader CR;
	CR.Open(InputFN);

	XProf XP;
	PDBChain Chain;
	uint Counter = 0;
	while (CR.GetNext(Chain))
		{
		if (++Counter%100 == 0)
			Progress("%s\r", Chain.m_Label.c_str());
		XP.Init(Chain);
		XP.ToCfv(g_fcfv);
		}
	Progress("\n");
	}

void XProf::Init(const PDBChain &Chain)
	{
	m_Chain = &Chain;
	m_L = m_Chain->GetSeqLength();
	}

void XProf::PosToCfv(FILE *f, uint Pos) const
	{
	if (f == 0)
		return;

	fprintf(f, "%u\t%c", Pos+1, m_Chain->m_Seq[Pos]);

	double x = m_Chain->m_Xs[Pos];
	double y = m_Chain->m_Ys[Pos];
	double z = m_Chain->m_Zs[Pos];

	fprintf(f, "\t%.1f\t%.1f\t%.1f", x, y, z);

	const uint FeatureCount = GetFeatureCount();
	for (uint FeatureIndex = 0; FeatureIndex < FeatureCount;
	  ++FeatureIndex)
		{
		double Value;
		uint iValue;
		GetFeature(FeatureIndex, Pos, Value, iValue);
		if (Value == DBL_MAX)
			fprintf(f, "\t.");
		else
			fprintf(f, "\t%.3g", Value);
		if (iValue == UINT_MAX)
			fprintf(f, "\t.");
		else
			fprintf(f, "\t%u", iValue);
		}
	fprintf(f, "\n");
	}

void XProf::ToCfv(FILE *f) const
	{
	if (f == 0)
		return;

	const uint L = GetSeqLength();
	const char *Label = m_Chain->m_Label.c_str();
	const uint FeatureCount = GetFeatureCount();

	fprintf(f, ">\t%s", Label);
	fprintf(f, "\t%u", L);
	fprintf(f, "\t%u", FeatureCount);
	for (uint i = 0; i < FeatureCount; ++i)
		{
		const char *Name = GetFeatureName(i);
		fprintf(f, "\t%s", Name);
		}
	fprintf(f, "\n");

	for (uint Pos = 0; Pos < L; ++Pos)
		PosToCfv(f, Pos);
	}

uint XProf::GetSphereNr(uint Pos, double Radius) const
	{
	vector<uint> PosVec;
	m_Chain->GetSphere(Pos, Radius, 0, m_L-1, PosVec);
	uint n = SIZE(PosVec);
	return n;
	}

void XProf::GetFeature(uint FeatureIndex, uint Pos,
  double &Value, uint &iValue) const
	{
	Value = DBL_MAX;
	iValue = UINT_MAX;
	switch (FeatureIndex)
		{
	case 0: Get_Ang_m2_p2(Pos, Value, iValue); return;
	case 1: Get_Ang_m3_p3(Pos, Value, iValue); return;
	case 2: Get_ED_p4(Pos, Value, iValue); return;
	case 3: Get_ED_m4(Pos, Value, iValue); return;
	case 4: Get_NU(Pos, Value, iValue); return;
	case 5: Get_ND(Pos, Value, iValue); return;
		}
	asserta(false);
	}

uint XProf::GetFeatureCount()
	{
	return 6;
	}

const char *XProf::GetFeatureName(uint FeatureIndex)
	{
	switch (FeatureIndex)
		{
	case 0: return "Ang_m2_p2";
	case 1: return "Ang_m3_p3";
	case 2: return "ED_p4";
	case 3: return "ED_m4";
	case 4: return "NU";
	case 5: return "ND";
		}
	asserta(false);
	return "???";
	}

double XProf::GetAngle(uint PosA1, uint PosA2, 
  uint PosB1, uint PosB2) const
	{
	vector<double> PtA1;
	vector<double> PtA2;
	vector<double> PtB1;
	vector<double> PtB2;

	m_Chain->GetPt(PosA1, PtA1);
	m_Chain->GetPt(PosA2, PtA2);

	m_Chain->GetPt(PosB1, PtB1);
	m_Chain->GetPt(PosB2, PtB2);

	vector<double> vA;
	vector<double> vB;
	Sub_Vecs(PtA2, PtA1, vA);
	Sub_Vecs(PtB2, PtB1, vB);

	double Radians = GetTheta_Vecs(vA, vB);
	double Deg = degrees(Radians);
	return Deg;
	}

void XProf::Get_NU(uint Pos, double &Value, uint &iValue) const
	{
	uint NU, ND;
	const uint L = GetSeqLength();
	GetHSE(*m_Chain, Pos, 12.0, 0, L-1, NU, ND);
	Value = double(NU);
	iValue = NU/3;
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_ND(uint Pos, double &Value, uint &iValue) const
	{
	uint NU, ND;
	const uint L = GetSeqLength();
	GetHSE(*m_Chain, Pos, 12.0, 0, L-1, NU, ND);
	Value = double(ND);
	iValue = ND;
	if (iValue > 10)
		iValue -= 10;
	iValue /= 2;
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang_m1_p1(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos == 0 || Pos + 1 >= m_L)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = GetAngle(Pos-1, Pos, Pos, Pos+1);
	iValue = uint(Value/20.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang_m2_p2(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos < 2 || Pos + 2 >= m_L)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = GetAngle(Pos-2, Pos, Pos, Pos+2);
	iValue = IntScale(Value, 1.0, 120.0);
	}

void XProf::Get_Ang_m3_p3(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos < 3 || Pos + 3 >= m_L)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = GetAngle(Pos-3, Pos, Pos, Pos+3);
	iValue = IntScale(Value, 1.0, 90.0);
	}

void XProf::Get_Ang01_23(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 3 >= m_L)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = GetAngle(Pos, Pos+1, Pos+2, Pos+3);
	iValue = uint(Value*4.0/80.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang01_34(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 4 >= m_L)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = GetAngle(Pos, Pos+1, Pos+3, Pos+4);
	iValue = uint(Value*4.0/80.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_Ang01_45(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 5 >= m_L)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = GetAngle(Pos, Pos+1, Pos+4, Pos+5);
	iValue = uint(Value*4.0/80.0);
	if (iValue >= XBINS)
		iValue = XBINS-1;
	}

void XProf::Get_ED_p4(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos + 4 >= m_L)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = m_Chain->GetDist(Pos, Pos+4);
	iValue = IntScale(Value, 4.5, 10.0);
	}

void XProf::Get_ED_m4(uint Pos, double &Value, uint &iValue) const
	{
	if (Pos < 4)
		{
		Value = DBL_MAX;
		iValue = UINT_MAX;
		return;
		}
	Value = m_Chain->GetDist(Pos-4, Pos);
	iValue = IntScale(Value, 4.5, 10.0);
	}

uint XProf::IntScale(double Value, double MinVal, double HiQ) const
	{
	double x = Value;
	if (x < MinVal)
		x = MinVal;
	x -= MinVal;
	uint i = (uint) ((x*XBINS)/(HiQ));
	if (i >= XBINS)
		i = XBINS - 1;
	return i;
	}
