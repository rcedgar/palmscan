#include "myutils.h"
#include "ppp.h"
#include "quarts.h"
#include "abcxyz.h"

double GetScoreNormal(double Mu, double Sigma, double x);

vector<string> PPP::m_FeatureNames;
uint PPP::m_ShapeLength = 64;
uint PPP::m_ShapeSmoothWindow = 2;
uint PPP::m_FeatureCount = 0;

const char *PPP::GetFeatureName(uint i)
	{
	asserta(i < SIZE(m_FeatureNames));
	return m_FeatureNames[i].c_str();
	}

void PPP::LogMe() const
	{
	Log("PPP::LogMe(), %u features, shape length %u, w %u, chains %u\n",
	  m_FeatureCount, m_ShapeLength, m_ShapeSmoothWindow, SIZE(m_FeatureValuesVec));

	for (uint i = 0; i < m_FeatureCount; ++i)
		{
		const char *Name = GetFeatureName(i);
		double Mean = m_FeatureMeans[i];
		double StdDev = m_FeatureStdDevs[i];

		Log("%6.6s", Name);

		if (Mean != DBL_MAX)
			Log("  Mean %.3g", Mean);

		if (Mean != DBL_MAX)
			Log("  StdDev %.3g", StdDev);

		Log("\n");
		}
	}

void PPP::ToPPP(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, ">%s\n", m_Name.c_str());
	for (uint i = 0; i < m_FeatureCount; ++i)
		{
		const char *Name = GetFeatureName(i);
		double Mean = m_FeatureMeans[i];
		double StdDev = m_FeatureStdDevs[i];
		fprintf(f, "%s\t%.4g\t%.4g\n",
		  Name, Mean, StdDev);
		}
	}

void PPP::FromPPP(const string &FileName)
	{
	Clear();

	FILE *f = OpenStdioFile(FileName);
	string Line;
	bool Ok = ReadLineStdioFile(f, Line);
	if (!Ok)
		Die("Empty PPP file '%s'", FileName.c_str());
	m_Name = Line.substr(1, string::npos);

	vector<string> Fields;
	for (uint FeatureIndex = 0; FeatureIndex < m_FeatureCount; ++FeatureIndex)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			Die("Truncated PPP file '%s'", FileName.c_str());

		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);

		const string &FeatureName = GetFeatureName(FeatureIndex);
		asserta(Fields[0] == FeatureName);
		double Mean = StrToFloat(Fields[1]);
		double StdDev = StrToFloat(Fields[2]);

		m_FeatureMeans.push_back(Mean);
		m_FeatureStdDevs.push_back(StdDev);
		}

	CloseStdioFile(f);
	}

void PPP::ToPDB(FILE *f) const
	{
	if (f == 0)
		return;

	fprintf(f, "TITLE     %s\n", m_Name.c_str());

#define C(x)	\
	fprintf(f, "REMARK    L" #x " mean %.3g, stddev %.3g\n",		\
		m_FeatureMeans[PF_L##x], m_FeatureStdDevs[PF_L##x]);	\
	fprintf(f, "REMARK    d" #x " mean %.3g, stddev %.3g\n",		\
		m_FeatureMeans[PF_d##x], m_FeatureStdDevs[PF_d##x]);

	C(AB)
	C(AC)
	C(BC)

#undef X

	for (uint i = 0; i < m_ShapeLength; ++i)
		{
		double Coord_X = m_FeatureMeans[PF_X00 + i];
		double Coord_Y = m_FeatureMeans[PF_Y00 + i];
		double Coord_Z = m_FeatureMeans[PF_Z00 + i];

		fprintf(f, "ATOM  ");			//  1 -  6        Record name   "ATOM  "
		fprintf(f, "%5u", i+1);			//  7 - 11        Integer       serial       Atom serial number.
		fprintf(f, " ");				// 12
		fprintf(f, " CA ");				// 13 - 16        Atom          name         Atom name.
		fprintf(f, " ");				// 17             Character     altLoc       Alternate location indicator.
		fprintf(f, "%3.3s", "ASP");		// 18 - 20        Residue name  resName      Residue name.
		fprintf(f, " ");				// 21
		fprintf(f, "%c", 'A');			// 22             Character     chainID      Chain identifier.
		fprintf(f, "%4u", i+1);			// 23 - 26        Integer       resSeq       Residue sequence number.
		fprintf(f, " ");				// 27             AChar         iCode        Code for insertion of residues.
		fprintf(f, "   ");				// 28 - 30
		fprintf(f, "%8.3f", Coord_X);	// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", Coord_Y);	// 39 - 46        Real(8.3)     y            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", Coord_Z);	// 47 - 57        Real(8.3)     z            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%6.2f", 1.0);		// 55 - 60        Real(6.2)     occupancy    Occupancy.
		fprintf(f, "%6.2f", 0.0);		// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		fprintf(f, "          ");		// 67 - 76
		fprintf(f, " C");				// 77 - 78        LString(2)    element      Element symbol, right-justified.
		fprintf(f, "  ");				// 79 - 80        LString(2)    charge       Charge on the atom.

		fprintf(f, "\n");
		}
	}

void PPP::Init(const string &Name)
	{
	Clear();
	m_Name = Name;
	m_FeatureMeans.resize(m_FeatureCount, DBL_MAX);
	m_FeatureStdDevs.resize(m_FeatureCount, DBL_MAX);
	}

void PPP::GlobalInit()
	{
	m_FeatureNames.clear();

#define F(x)	m_FeatureNames.push_back(#x);
#include "pppfeatures.h"

	m_FeatureCount = SIZE(m_FeatureNames);
	}

void PPP::AddChain(const PDBChain &Chain)
	{
	vector<double> FeatureValues;
	GetFeatureVec(Chain, FeatureValues);
	asserta(SIZE(FeatureValues) == m_FeatureCount);
	m_FeatureValuesVec.push_back(FeatureValues);
	}

void PPP::GetFeatureVec(const PDBChain &Chain,
  vector<double> &FeatureValues) const
	{
	FeatureValues.clear();
	FeatureValues.resize(m_FeatureCount, DBL_MAX);

	Chain.CheckPPCMotifCoords();

	uint PosA = Chain.GetMotifPos(A);
	uint PosB = Chain.GetMotifPos(B);
	uint PosC = Chain.GetMotifPos(C);

	FeatureValues[PF_LAB] = PosB - PosA + 1;
	FeatureValues[PF_LAC] = PosC - PosA + 1;
	FeatureValues[PF_LBC] = PosC - PosB + 1;

	FeatureValues[PF_dAB] = Chain.GetDist(PosA, PosB);
	FeatureValues[PF_dAC] = Chain.GetDist(PosA, PosC);
	FeatureValues[PF_dBC] = Chain.GetDist(PosB, PosC);

	for (uint i = 0; i < m_ShapeLength; ++i)
		{
		double Coord_X = 
		  Chain.GetSmoothedCoord(X, i, m_ShapeLength, m_ShapeSmoothWindow);

		double Coord_Y = 
		  Chain.GetSmoothedCoord(Y, i, m_ShapeLength, m_ShapeSmoothWindow);

		double Coord_Z = 
		  Chain.GetSmoothedCoord(Z, i, m_ShapeLength, m_ShapeSmoothWindow);

		asserta(PF_X00 + i < m_FeatureCount);
		asserta(PF_Y00 + i < m_FeatureCount);
		asserta(PF_Z00 + i < m_FeatureCount);

		FeatureValues[PF_X00 + i] = Coord_X;
		FeatureValues[PF_Y00 + i] = Coord_Y;
		FeatureValues[PF_Z00 + i] = Coord_Z;
		}
	}

void PPP::Averages()
	{
	for (uint i = 0; i < m_FeatureCount; ++i)
		Average(i);
	}

void PPP::Average(uint FeatureIndex)
	{
	const uint ChainCount = SIZE(m_FeatureValuesVec);
	asserta(ChainCount > 0);
	vector<double> Values;
	for (uint ChainIndex = 0; ChainIndex< ChainCount; ++ChainIndex)
		{
		asserta(ChainIndex < SIZE(m_FeatureValuesVec));
		asserta(FeatureIndex < SIZE(m_FeatureValuesVec[ChainIndex]));
		double Value = m_FeatureValuesVec[ChainIndex][FeatureIndex];
		Values.push_back(Value);
		}
	QuartsDouble QD;
	GetQuartsDouble(Values, QD);
	asserta(FeatureIndex < SIZE(m_FeatureMeans));
	asserta(FeatureIndex < SIZE(m_FeatureStdDevs));

	m_FeatureMeans[FeatureIndex] = QD.Avg;
	m_FeatureStdDevs[FeatureIndex] = QD.StdDev;
	}

double PPP::ScoreChain(const PDBChain &Chain,
  uint PosA, uint PosB, uint PosC) const
	{
	vector<vector<double> > MotifCoords(3);

	Chain.GetPt(PosA, MotifCoords[A]);
	Chain.GetPt(PosB, MotifCoords[B]);
	Chain.GetPt(PosC, MotifCoords[C]);

	vector<double> t;
	vector<vector<double> > R;
	GetTriForm(MotifCoords, t, R);

	vector<double> Pt(3);
	vector<double> XPt;
	double SumScore = 0;
	for (uint i = 0; i < m_ShapeLength; ++i)
		{
		Pt[X] = Chain.GetSmoothedCoord(X, i, m_ShapeLength, m_ShapeSmoothWindow);
		Pt[Y] = Chain.GetSmoothedCoord(Y, i, m_ShapeLength, m_ShapeSmoothWindow);
		Pt[Z] = Chain.GetSmoothedCoord(Z, i, m_ShapeLength, m_ShapeSmoothWindow);
		XFormPt(Pt, t, R, XPt);

		double ChainX = XPt[X];
		double ChainY = XPt[Y];
		double ChainZ = XPt[Z];

		double MeanX = m_FeatureMeans[PF_X00 + i];
		double MeanY = m_FeatureMeans[PF_Y00 + i];
		double MeanZ = m_FeatureMeans[PF_Z00 + i];

		double StdDevX = m_FeatureStdDevs[PF_X00 + i];
		double StdDevY = m_FeatureStdDevs[PF_Y00 + i];
		double StdDevZ = m_FeatureStdDevs[PF_Z00 + i];

		double ScoreX = GetScoreNormal(MeanX, StdDevX, ChainX);
		double ScoreY = GetScoreNormal(MeanY, StdDevY, ChainY);
		double ScoreZ = GetScoreNormal(MeanZ, StdDevZ, ChainZ);
//		SumScore += ScoreX*ScoreX + ScoreY*ScoreY + ScoreZ*ScoreZ;
		SumScore += ScoreX + ScoreY + ScoreZ;
		}

	double Score = SumScore/(m_ShapeLength*3);
	return Score;
	}
