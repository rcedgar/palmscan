#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

void GetThreeFromOne(char aa, string &AAA);

void ReadLinesFromFile(const string &FileName, vector<string> &Lines)
	{
	Lines.clear();
	FILE *f = OpenStdioFile(FileName);
	string Line;
	while (ReadLineStdioFile(f, Line))
		{
		Lines.push_back(Line);
		}
	CloseStdioFile(f);
	}

/***
PDBChain ATOM record format
http://www.wwpdb.org/documentation/file-format-content/format33/sect9.html#ATOM

COLUMNS        DATA  TYPE    FIELD        DEFINITION
-------------------------------------------------------------------------------------
 1 -  6        Record name   "ATOM  "
 7 - 11        Integer       serial       Atom  serial number.
13 - 16        Atom          name         Atom name.
17             Character     altLoc       Alternate location indicator.
18 - 20        Residue name  resName      Residue name.
22             Character     chainID      Chain identifier.
23 - 26        Integer       resSeq       Residue sequence number.
27             AChar         iCode        Code for insertion of residues.
31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
47 - 54        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
55 - 60        Real(6.2)     occupancy    Occupancy.
61 - 66        Real(6.2)     tempFactor   Temperature  factor.
77 - 78        LString(2)    element      Element symbol, right-justified.
79 - 80        LString(2)    charge       Charge  on the atom.

MODEL     1                                                                     

          1         2         3         4         5         6         7         
01234567890123456789012345678901234567890123456789012345678901234567890123456789
                              xxxxxxxxyyyyyyyyzzzzzzzz
ATOM      1  N   PHE A   1      34.582  19.022  -8.646  1.00 24.35           N  
ATOM      2  CA  PHE A   1      33.319  19.558  -8.153  1.00 24.35           C  
ATOM      3  C   PHE A   1      32.243  19.483  -9.229  1.00 24.35           C  
ATOM      4  CB  PHE A   1      33.494  21.006  -7.685  1.00 24.35           C  
***/

void PDBChain::LogMe(bool WithCoords) const
	{
	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	Log("\n");
	Log(">%s\n", m_Label.c_str());
	Log("%s\n", m_Seq.c_str());
	if (!m_MotifPosVec.empty())
		{
		asserta(m_MotifPosVec.size() == 3);
		string MotifA, MotifB, MotifC;
		GetMotifSeq(A, MotifA);
		GetMotifSeq(B, MotifB);
		GetMotifSeq(C, MotifC);

		Log("A:%u:%s  B:%u:%s  C:%u:%s\n",
		  m_MotifPosVec[0]+1, MotifA.c_str(),
		  m_MotifPosVec[1]+1, MotifB.c_str(),
		  m_MotifPosVec[2]+1, MotifC.c_str());
		}

	if (WithCoords)
		{
		Log("\n");
		Log(" Pos  S         X         Y         Z\n");
		for (size_t i = 0; i < L; ++i)
			{
			Log("%4u", i);
			Log("  %c", m_Seq[i]);
			Log("  %8.3f", m_Xs[i]);
			Log("  %8.3f", m_Ys[i]);
			Log("  %8.3f", m_Zs[i]);
			Log("\n");
			}
		}
	}

void PDBChain::ToCal(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	ToCal(f);
	CloseStdioFile(f);
	}

void PDBChain::ToCal(FILE* f) const
	{
	if (f == 0)
		return;
	uint QL = GetSeqLength();
	ToCalSeg(f, 0, QL);
	}

void PDBChain::ToCalSeg(FILE *f, uint Pos, uint n) const
	{
	if (f == 0)
		return;
	if (n == 0)
		return;

	const size_t L = m_Xs.size();
	if (L == 0)
		return;

	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	if (m_MotifPosVec.empty())
		fprintf(f, ">%s\n", m_Label.c_str());
	else
		{
		asserta(m_MotifPosVec.size() == 3);

		vector<string> Fields;
		Split(m_Label, Fields, ' ');
		string NewLabel = Fields[0];

		uint PosA = m_MotifPosVec[0];
		uint PosB = m_MotifPosVec[1];
		uint PosC = m_MotifPosVec[2];
		asserta(PosA >= Pos);
		asserta(PosC + CL <= Pos + n);

		string MotifA, MotifB, MotifC;
		GetMotifSeq(A, MotifA);
		GetMotifSeq(B, MotifB);
		GetMotifSeq(C, MotifC);

		fprintf(f, ">%s A:%u:%s B:%u:%s C:%u:%s\n",
		  NewLabel.c_str(),
		  m_MotifPosVec[0]+1, MotifA.c_str(),
		  m_MotifPosVec[1]+1, MotifB.c_str(),
		  m_MotifPosVec[2]+1, MotifC.c_str());
		}

	for (size_t i = Pos; i < Pos + n; ++i)
		{
		if (i < SIZE(m_Seq))
			fputc(m_Seq[i], f);
		else
			fputc('*', f);
		fprintf(f, "\t%.3f", m_Xs[i]);
		fprintf(f, "\t%.3f", m_Ys[i]);
		fprintf(f, "\t%.3f", m_Zs[i]);
		fputc('\n', f);
		}
	}

void PDBChain::ToPDB(const string &FileName) const
	{
	if (FileName == "")
		return;

	FILE *f = CreateStdioFile(FileName);
	fprintf(f, "TITLE %s\n", m_Label.c_str());

	uint AtomCount = SIZE(m_ATOMs);
	if (AtomCount > 0)
		{
		for (uint i = 0; i < AtomCount; ++i)
			{
			const vector<string> &v = m_ATOMs[i];
			for (uint j = 0; j < SIZE(v); ++j)
				{
				fputs(v[j].c_str(), f);
				fputc('\n', f);
				}
			}
		return;
		}

	const size_t L = m_Xs.size();
	asserta(m_Ys.size() == L);
	asserta(m_Zs.size() == L);
	asserta(m_Seq.size() == L);

	const char Chain = 'A';

	for (uint i = 0; i < L; ++i)
		{
		char aa = m_Seq[i];
		string sAAA;
		GetThreeFromOne(aa, sAAA);
		const char *AAA = sAAA.c_str();

		fprintf(f, "ATOM  ");			//  1 -  6        Record name   "ATOM  "
		fprintf(f, "%5u", i+1);			//  7 - 11        Integer       serial       Atom serial number.
		fprintf(f, " ");				// 12
		fprintf(f, " CA ");				// 13 - 16        Atom          name         Atom name.
		fprintf(f, " ");				// 17             Character     altLoc       Alternate location indicator.
		fprintf(f, "%3.3s", AAA);		// 18 - 20        Residue name  resName      Residue name.
		fprintf(f, " ");				// 21
		fprintf(f, "%c", Chain);		// 22             Character     chainID      Chain identifier.
		fprintf(f, "%4u", i+1);			// 23 - 26        Integer       resSeq       Residue sequence number.
		fprintf(f, " ");				// 27             AChar         iCode        Code for insertion of residues.
		fprintf(f, "   ");				// 28 - 30
		fprintf(f, "%8.3f", m_Xs[i]);	// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
		fprintf(f, "%8.3f", m_Ys[i]);	// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
		fprintf(f, "%8.3f", m_Zs[i]);	// 47 - 57        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
		fprintf(f, "%6.2f", 1.0);		// 55 - 60        Real(6.2)     occupancy    Occupancy.
		fprintf(f, "%6.2f", 0.0);		// 61 - 66        Real(6.2)     tempFactor   Temperature  factor.
		fprintf(f, "          ");		// 67 - 76
		fprintf(f, " C");				// 77 - 78        LString(2)    element      Element symbol, right-justified.
		fprintf(f, "  ");				// 79 - 80        LString(2)    charge       Charge on the atom.

		fprintf(f, "\n");
		}
	}

void PDBChain::RenumberResidues(uint Start)
	{
	uint N = SIZE(m_ATOMs);
	asserta(N > 0);
	uint CurrInputResidueNr = UINT_MAX;
	vector<vector<string> > OutputAtoms;
	for (uint OutputResidueIndex = 0; OutputResidueIndex < N;
	  ++OutputResidueIndex)
		{
		const vector<string> &InputLines = m_ATOMs[OutputResidueIndex];
		vector<string> OutputLines;
		const uint M = SIZE(InputLines);
		asserta(M > 0);
		for (uint j = 0; j < M; ++j)
			{
			const string &InputLine = InputLines[j];
			uint InputResidueNr = GetResidueNrFromATOMLine(InputLine);
			string OutputLine;
			SetResidueNrInATOMLine(InputLine,
			  Start+OutputResidueIndex+1, OutputLine);
			OutputLines.push_back(OutputLine);
			}
		OutputAtoms.push_back(OutputLines);
		}
	m_ATOMs = OutputAtoms;
	}

static double GetFloatFromString(const string &s, uint Pos, uint n)
	{
	string t = s.substr(Pos, n);
	StripWhiteSpace(t);
	double Value = StrToFloat(t);
	return Value;
	}

// 31 - 38        Real(8.3)     x            Orthogonal coordinates for X in Angstroms.
// 39 - 46        Real(8.3)     y            Orthogonal coordinates for Y in Angstroms.
// 47 - 57        Real(8.3)     z            Orthogonal coordinates for Z in Angstroms.
void PDBChain::GetXYZFromATOMLine(const string &InputLine,
  double &x, double &y, double &z)
	{
	x = GetFloatFromString(InputLine, 30, 8);
	y = GetFloatFromString(InputLine, 38, 8);
	z = GetFloatFromString(InputLine, 46, 8);
	}

void PDBChain::SetXYZInATOMLine(const string &InputLine,
  double x, double y, double z, string &OutputLine)
	{
	string sx;
	string sy;
	string sz;
	Ps(sx, "%8.3f", x);
	Ps(sy, "%8.3f", y);
	Ps(sz, "%8.3f", z);
	asserta(SIZE(sx) == 8);
	asserta(SIZE(sy) == 8);
	asserta(SIZE(sz) == 8);

	OutputLine = InputLine;
	for (uint i = 0; i < 8; ++i)
		{
		OutputLine[30+i] = sx[i];
		OutputLine[38+i] = sy[i];
		OutputLine[46+i] = sz[i];
		}
	}

void PDBChain::GetCAAtomLine(uint Pos, string &Line) const
	{
	asserta(Pos < SIZE(m_ATOMs));
	const vector<string> &v = m_ATOMs[Pos];
	for (uint i = 0; i < SIZE(v); ++i)
		{
		Line = v[i];
		//asserta(SIZE(Line) > 40);
		//string AtomName = Line.substr(12, 4);
		//StripWhiteSpace(AtomName);
		string AtomName;
		GetAtomNameFromATOMLine(Line, AtomName);
		if (AtomName == "CA")
			return;
		}
	assert(false);
	}

int PDBChain::GetResidueNr(uint Pos) const
	{
	//string CALine;
	//GetCAAtomLine(Pos, CALine);
	//int ResNr = GetResidueNrFromATOMLine(CALine);
	if (SIZE(m_ResNrs) == 0)
		return INT_MAX;
	asserta(Pos < SIZE(m_ResNrs));
	return m_ResNrs[Pos];
	}

void PDBChain::GetResidueRange(uint PosLo, uint ResidueCount,
  int &ResLo, int &ResHi) const
	{
	ResLo = INT_MAX;
	ResHi = INT_MIN;
	for (uint Pos = PosLo; Pos < PosLo + ResidueCount; ++Pos)
		{
		int ResNr = GetResidueNr(Pos);
		if (ResNr < ResLo)
			ResLo = ResNr;
		if (ResNr > ResHi)
			ResHi = ResNr;
		}
	}

void PDBChain::GetAtomNameFromATOMLine(const string &Line,
  string &AtomName)
	{
	asserta(SIZE(Line) > 40);
	AtomName = Line.substr(12, 4);
	StripWhiteSpace(AtomName);
	}

// Residue nr in Cols 23-26 (1-based)
int PDBChain::GetResidueNrFromATOMLine(const string &Line)
	{
	asserta(SIZE(Line) > 60);
	string s;
	s += Line[22];
	s += Line[23];
	s += Line[24];
	s += Line[25];
	StripWhiteSpace(s);
	int Nr = atoi(s.c_str());
	return Nr;
	}

// Residue nr in Cols 23-26 (1-based)
void PDBChain::SetResidueNrInATOMLine(const string &InputLine,
  uint ResidueNr, string &OutputLine)
	{
	string s;
	Ps(s, "%4u", ResidueNr);
	uint n = SIZE(s);
	if (n != 4)
		Die("Residue number %d overflow (max 9999)", ResidueNr);
	OutputLine.clear();
	uint N = SIZE(InputLine);
	for (uint i = 0; i < 23; ++i)
		OutputLine += InputLine[i];
	for (uint i = 0; i < 4; ++i)
		OutputLine += s[i];
	for (uint i = 23+4; i < N; ++i)
		OutputLine += InputLine[i];
	asserta(SIZE(OutputLine) == SIZE(InputLine));
	}

void PDBChain::HackHETAMLine(string &Line)
	{
	if (strncmp(Line.c_str(), "HETATM", 6) == 0 &&
		Line.substr(17, 3) == "MSE")
		{
		Line[0] = 'A';
		Line[1] = 'T';
		Line[2] = 'O';
		Line[3] = 'M';
		Line[4] = ' ';
		Line[5] = ' ';

		Line[17] = 'M';
		Line[18] = 'E';
		Line[19] = 'T';
		}
	}

char PDBChain::FromPDBLines(const string &Label,
  const vector<string> &Lines)
	{
	Clear();

	m_Label = Label;
	const uint N = SIZE(Lines);
	char Chain = 0;
	uint ResidueCount = 0;
	int CurrentResidueNumber = INT_MAX;
	for (uint LineNr = 0; LineNr < N; ++LineNr)
		{
		string Line = Lines[LineNr];
		const size_t L = Line.size();

		HackHETAMLine(Line);
		if (strncmp(Line.c_str(), "ATOM  ", 6) != 0)
			continue;
		if (Line[26] != ' ') // Insertion code
			continue;

		char LineChain = Line[21];
		if (Chain == 0)
			Chain = LineChain;
		else if (Chain != LineChain)
			Die("PDBChain::FromLines() two chains %c, %c",
			  Chain, LineChain);

		int ResidueNumber = GetResidueNrFromATOMLine(Line);
		if (ResidueNumber != CurrentResidueNumber)
			{
			CurrentResidueNumber = ResidueNumber;
			++ResidueCount;
			m_ATOMs.resize(ResidueCount);
			}
		asserta(SIZE(m_ATOMs) > 0);
		m_ATOMs.back().push_back(Line);
		}

	AppendChainToLabel(m_Label, Chain);

	const uint L = SIZE(m_ATOMs);
	for (uint i = 0; i < L; ++i)
		{
		char aa;
		double X, Y, Z;
		int ResNr;
		bool CAFound = 
		  GetFieldsFromResidueATOMLines(m_ATOMs[i], X, Y, Z, aa, ResNr);
		if (!CAFound)
			continue;

		m_Seq.push_back(aa);
		m_Xs.push_back(X);
		m_Ys.push_back(Y);
		m_Zs.push_back(Z);
		m_ResNrs.push_back(ResNr);
		}

	return Chain;
	}

bool PDBChain::GetFieldsFromResidueATOMLines(const vector<string> &Lines,
  double &X, double &Y, double &Z, char &aa, int &ResNr)
	{
	aa = 'X';
	ResNr = -999;
	X = -999;
	Y = -999;
	Z = -999;
	const uint N = SIZE(Lines);
	for (uint i = 0; i < SIZE(Lines); ++i)
		{
		const string &Line = Lines[i];
		string AtomName = Line.substr(12, 4);
		StripWhiteSpace(AtomName);
		if (AtomName == "CA")
			{
			string AAA = Line.substr(17, 3);
			aa = GetOneFromThree(AAA);

			string sX, sY, sZ;
			sX = Line.substr(30, 8);
			sY = Line.substr(38, 8);
			sZ = Line.substr(46, 8);

			StripWhiteSpace(sX);
			StripWhiteSpace(sY);
			StripWhiteSpace(sZ);

			X = StrToFloat(sX);
			Y = StrToFloat(sY);
			Z = StrToFloat(sZ);

			ResNr = GetResidueNrFromATOMLine(Line);
			return true;
			}
		}
	return false;
	}

void PDBChain::GetSubSeq(uint Pos, uint n, string &s) const
	{
	s.clear();
	size_t L = m_Seq.size();
	asserta(Pos + n <= L);

	for (uint i = 0; i < n; ++i)
		{
		char c = m_Seq[Pos+i];
		s += c;
		}
	}

uint PDBChain::GetMotifLength(uint i)
	{
	asserta(i < 3);
	static const uint Lengths[3] = { AL, BL, CL };
	return Lengths[i];
	}

void PDBChain::GetMotifSeq(uint i, string &Seq) const
	{
	asserta(m_MotifPosVec.size() == 3);
	uint Pos = GetMotifPos(i);
	uint ML = GetMotifLength(i);
	GetSubSeq(Pos, ML, Seq);
	}

double PDBChain::GetCoord(uint Axis, uint Pos) const
	{
	switch (Axis)
		{
	case X: return m_Xs[Pos];
	case Y: return m_Ys[Pos];
	case Z: return m_Zs[Pos];
		}
	asserta(false);
	return 0;
	}

void PDBChain::GetPt(uint Pos, vector<double> &Pt) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	Resize3(Pt);
	Pt[X] = m_Xs[Pos];
	Pt[Y] = m_Ys[Pos];
	Pt[Z] = m_Zs[Pos];
	}

void PDBChain::SetPt(uint Pos, const vector<double> &Pt)
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));

	m_Xs[Pos] = Pt[X];
	m_Ys[Pos] = Pt[Y];
	m_Zs[Pos] = Pt[Z];
	}

bool PDBChain::Get_CB_Pt(uint Pos, vector<double> &Pt) const
	{
	Pt.resize(3);
	bool Ok = Get_CB_XYZ(Pos, Pt[0], Pt[1], Pt[2]);
	return Ok;
	}

bool PDBChain::Get_CB_XYZ(uint Pos, double &x, double &y, double &z) const
	{
	x = DBL_MAX;
	y = DBL_MAX;
	z = DBL_MAX;
	if (m_ATOMs.empty())
		return false;
	asserta(Pos < SIZE(m_ATOMs));
	const vector<string> &v = m_ATOMs[Pos];
	for (uint i = 0; i < SIZE(v); ++i)
		{
		const string &Line = v[i];
		string AtomName;
		GetAtomNameFromATOMLine(Line, AtomName);
		if (AtomName == "CB")
			{
			GetXYZFromATOMLine(Line, x, y, z);
			return true;
			}
		}
	return false;
	}

void PDBChain::GetXYZ(uint Pos, double &x, double &y, double &z) const
	{
	assert(Pos < SIZE(m_Xs));
	assert(Pos < SIZE(m_Ys));
	assert(Pos < SIZE(m_Zs));
	x = m_Xs[Pos];
	y = m_Ys[Pos];
	z = m_Zs[Pos];
	}

void PDBChain::GetDistMx(uint Pos, uint L,
  vector<vector<double> > &Mx) const
	{
	Mx.resize(L);
	for (uint i = 0; i < L; ++i)
		{
		Mx[i].resize(L);
		for (uint j = 0; j < L; ++j)
			{
			double d = GetDist(Pos+i, Pos+j);
			Mx[i][j] = d;
			}
		}
	}

double PDBChain::GetDist(uint Pos1, uint Pos2) const
	{
	double x1, y1, z1;
	double x2, y2, z2;
	GetXYZ(Pos1, x1, y1, z1);
	GetXYZ(Pos2, x2, y2, z2);
	double d = GetDist3D(x1, y1, z1, x2, y2, z2);
	return d;
	}

double PDBChain::GetDist2(uint Pos1, uint Pos2) const
	{
	double x1 = m_Xs[Pos1];
	double y1 = m_Ys[Pos1];
	double z1 = m_Zs[Pos1];

	double x2 = m_Xs[Pos2];
	double y2 = m_Ys[Pos2];
	double z2 = m_Zs[Pos2];

	double dx = x1 - x2;
	double dy = y1 - y2;
	double dz = z1 - z2;

	double d = GetDist(Pos1, Pos2);

	double d2 = dx*dx + dy*dy + dz*dz;
	asserta(feq(d*d, d2));
	return d2;
	}

void PDBChain::GetMotifDists(double &AB, double &BC, double &AC) const
	{
	vector<vector<double> > MotifCoords;
	GetMotifCoords(MotifCoords);
	AB = GetDist_Mxij(MotifCoords, A, B);
	BC = GetDist_Mxij(MotifCoords, B, C);
	AC = GetDist_Mxij(MotifCoords, A, C);
	}

void PDBChain::GetMotifDists2(double &AB, double &BC, double &AC) const
	{
	vector<vector<double> > MotifCoords;
	GetMotifCoords(MotifCoords);
	AB = GetDist2_Mxij(MotifCoords, A, B);
	BC = GetDist2_Mxij(MotifCoords, B, C);
	AC = GetDist2_Mxij(MotifCoords, A, C);
	}

void PDBChain::GetDGDCoords(vector<vector<double> > &DGDCoords) const
	{
	asserta(m_MotifPosVec.size() == 3);

	Resize3x3(DGDCoords);

	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	uint L = GetSeqLength();
	uint PosA_D = PosA + 3;
	uint PosB_G = PosB + 1;
	uint PosC_D = PosC + 3;
	asserta(PosA_D < PosB_G);
	asserta(PosB_G < PosC_D);
	asserta(PosC_D < L);

	asserta(m_Seq[PosA_D] == 'D');
	asserta(m_Seq[PosB_G] == 'G');
	asserta(m_Seq[PosC_D] == 'D');

	DGDCoords[A][X] = m_Xs[PosA_D];
	DGDCoords[A][Y] = m_Ys[PosA_D];
	DGDCoords[A][Z] = m_Zs[PosA_D];

	DGDCoords[B][X] = m_Xs[PosB_G];
	DGDCoords[B][Y] = m_Ys[PosB_G];
	DGDCoords[B][Z] = m_Zs[PosB_G];

	DGDCoords[C][X] = m_Xs[PosC_D];
	DGDCoords[C][Y] = m_Ys[PosC_D];
	DGDCoords[C][Z] = m_Zs[PosC_D];
	}

void PDBChain::GetMotifCoords(vector<vector<double> > &MotifCoords) const
	{
	asserta(m_MotifPosVec.size() == 3);

	Resize3x3(MotifCoords);

	uint PosA = m_MotifPosVec[0];
	uint PosB = m_MotifPosVec[1];
	uint PosC = m_MotifPosVec[2];

	MotifCoords[A][X] = m_Xs[PosA];
	MotifCoords[A][Y] = m_Ys[PosA];
	MotifCoords[A][Z] = m_Zs[PosA];

	MotifCoords[B][X] = m_Xs[PosB];
	MotifCoords[B][Y] = m_Ys[PosB];
	MotifCoords[B][Z] = m_Zs[PosB];

	MotifCoords[C][X] = m_Xs[PosC];
	MotifCoords[C][Y] = m_Ys[PosC];
	MotifCoords[C][Z] = m_Zs[PosC];
	}

void PDBChain::SetMotifPosVec(uint PosA, uint PosB, uint PosC)
	{
	const uint L = GetSeqLength();
	asserta(PosA < L);
	asserta(PosB < L);
	asserta(PosC < L);
	asserta(PosA < PosB);
	asserta(PosC < PosA || PosB < PosC);
	m_MotifPosVec.clear();
	m_MotifPosVec.push_back(PosA);
	m_MotifPosVec.push_back(PosB);
	m_MotifPosVec.push_back(PosC);
	}

void PDBChain::PrintSeqCoords(FILE *f) const
	{
	if (f == 0)
		return;
	const uint L = GetSeqLength();

	fprintf(f, "\n");
	fprintf(f, ">%s\n", m_Label.c_str());
	for (uint i = 0; i < L; ++i)
		{
		if ((i+1)%10 == 0)
			fprintf(f, "%10u", (i+1)/10);
		}
	fprintf(f, "\n");

	for (uint i = 0; i < L; ++i)
		fprintf(f, "%u", (i+1)%10);
	fprintf(f, "\n");

	fprintf(f, "%s\n", m_Seq.c_str());
	}

uint PDBChain::GetSeqLength() const
	{
	return SIZE(m_Seq);
	}

uint PDBChain::GetMotifPos(uint MotifIndex) const
	{
	asserta(MotifIndex < 3);
	asserta(m_MotifPosVec.size() == 3);
	uint Pos = m_MotifPosVec[MotifIndex];
	return Pos;
	}

void PDBChain::GetSubSeq(uint MotifStartPos, uint n,
  bool FailOnOverflow, string &MotifSeq) const
	{
	MotifSeq.clear();
	const uint L = GetSeqLength();
	for (uint i = 0; i < n; ++i)
		{
		uint Pos = MotifStartPos + i;
		if (Pos < 0 || Pos >= L)
			{
			if (FailOnOverflow)
				Die("'%s' GetMotifSeqFromMidPos overflow",
				  m_Label.c_str());
			MotifSeq += '!';
			}
		else
			MotifSeq += m_Seq[Pos];
		}
	}

void PDBChain::ChainsFromLines(const string &Label,
  const vector<string> &Lines, vector<PDBChain *> &Chains)
	{
	Chains.clear();
	const uint N = SIZE(Lines);
	vector<string> ChainLines;
	char CurrChainChar = 0;
	bool AnyAtoms = false;
	for (uint i = 0; i < N; ++i)
		{
		const string &Line = Lines[i];

		if (IsATOMLine(Line))
			{
			if (Line.size() < 57)
				continue;
			char ChainChar = Line[21];
			if (ChainChar != CurrChainChar)
				{
				if (AnyAtoms && !ChainLines.empty())
					{
					PDBChain *Chain = new PDBChain;
					char ChainChar = Chain->FromPDBLines(Label, ChainLines);
					if (ChainChar != 0)
						Chains.push_back(Chain);
					ChainLines.clear();
					AnyAtoms = false;
					}
				CurrChainChar = ChainChar;
				}
			ChainLines.push_back(Line);
			AnyAtoms = true;
			}
		}

	if (!ChainLines.empty() && AnyAtoms)
		{
		PDBChain *Chain = new PDBChain;
		Chain->FromPDBLines(Label, ChainLines);
		ChainLines.clear();
		Chains.push_back(Chain);
		}
	}

void PDBChain::ReadChainsFromFile(const string &FileName,
  vector<PDBChain *> &Chains)
	{
	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);
	ChainsFromLines(FileName, Lines, Chains);
	}

void PDBChain::GetTriFormChain_DGD(PDBChain &XChain) const
	{
	vector<vector<double> > DGDCoords;
	GetDGDCoords(DGDCoords);

	vector<double> t;
	vector<vector<double> > R;
	GetTriForm(DGDCoords, t, R);
	GetTriFormChain_tR(t, R, XChain);

#if 1
	{
	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	uint PosA_D = PosA + 3;
	uint PosB_G = PosB + 1;
	uint PosC_D = PosC + 3;

	Log("\n");
	Log("GetTriFormChain_DGD()\n");

	Log("  In  MotifA_D");
	Log("  %c[%3u]", m_Seq[PosA_D], PosA_D);
	Log("  %8.3f", m_Xs[PosA_D]);
	Log("  %8.3f", m_Ys[PosA_D]);
	Log("  %8.3f", m_Zs[PosA_D]);
	Log("\n");

	Log("  In  MotifB_G");
	Log("  %c[%3u]", m_Seq[PosB_G], PosB_G);
	Log("  %8.3f", m_Xs[PosB_G]);
	Log("  %8.3f", m_Ys[PosB_G]);
	Log("  %8.3f", m_Zs[PosB_G]);
	Log("\n");

	Log("  In  MotifC_D");
	Log("  %c[%3u]", m_Seq[PosC_D], PosC_D);
	Log("  %8.3f", m_Xs[PosC_D]);
	Log("  %8.3f", m_Ys[PosC_D]);
	Log("  %8.3f", m_Zs[PosC_D]);
	Log("\n");

	Log("\n");
	Log(" Out  MotifA_D");
	Log("  %c[%3u]", XChain.m_Seq[PosA_D], PosA_D);
	Log("  %8.3f", XChain.m_Xs[PosA_D]);
	Log("  %8.3f", XChain.m_Ys[PosA_D]);
	Log("  %8.3f", XChain.m_Zs[PosA_D]);
	Log("\n");

	Log(" Out  MotifB_G");
	Log("  %c[%3u]", XChain.m_Seq[PosB_G], PosB_G);
	Log("  %8.3f", XChain.m_Xs[PosB_G]);
	Log("  %8.3f", XChain.m_Ys[PosB_G]);
	Log("  %8.3f", XChain.m_Zs[PosB_G]);
	Log("\n");

	Log(" Out  MotifC_D");
	Log("  %c[%3u]", XChain.m_Seq[PosC_D], PosC_D);
	Log("  %8.3f", XChain.m_Xs[PosC_D]);
	Log("  %8.3f", XChain.m_Ys[PosC_D]);
	Log("  %8.3f", XChain.m_Zs[PosC_D]);
	Log("\n");
	}
#endif
	}

void PDBChain::GetTriFormChain_MotifCoords(PDBChain &XChain) const
	{
	vector<vector<double> > MotifCoords;
	GetMotifCoords(MotifCoords);

	vector<double> t;
	vector<vector<double> > R;
	GetTriForm(MotifCoords, t, R);

	GetTriFormChain_tR(t, R, XChain);
	}

void PDBChain::GetTriFormChain_tR(
  const vector<double> &t,
  const vector<vector<double> > &R,
  PDBChain &XChain) const
	{
	XChain.Clear();
	XChain.m_Label = m_Label;
	XChain.m_Seq = m_Seq;
	XChain.m_MotifPosVec = m_MotifPosVec;
	XChain.m_ATOMs = m_ATOMs;

	const uint N = SIZE(m_Seq);
	asserta(SIZE(m_Xs) == N);
	asserta(SIZE(m_Ys) == N);
	asserta(SIZE(m_Zs) == N);

	const uint AtomCount = SIZE(m_ATOMs);
	asserta(AtomCount == 0 || AtomCount == N);

	vector<double> Pt(3);
	vector<double> XPt(3);
	for (uint Pos = 0; Pos < N; ++Pos)
		{
		GetPt(Pos, Pt);
		XFormPt(Pt, t, R, XPt);

		double x = XPt[X];
		double y = XPt[Y];
		double z = XPt[Z];

		XChain.m_Xs.push_back(x);
		XChain.m_Ys.push_back(y);
		XChain.m_Zs.push_back(z);

		if (AtomCount != 0)
			{
			vector<string> &a = XChain.m_ATOMs[Pos];
			for (uint i = 0; i < SIZE(a); ++i)
				{
				GetXYZFromATOMLine(a[i], Pt[X], Pt[Y], Pt[Z]);
				XFormPt(Pt, t, R, XPt);
				SetXYZInATOMLine(a[i], XPt[X], XPt[Y], XPt[Z], a[i]);
				}
			}
		}
	}

void PDBChain::AppendChainToLabel(string &Label, char Chain)
	{
	if (Chain == 0)
		return;

	string _X = "_";
	_X += Chain;
	if (Label.find(_X) == string::npos)
		{
		Label += "_";
		Label += Chain;
		}
	}

uint PDBChain::GetPalmPrintLength(uint PosA, uint PosC, uint L)
	{
	asserta(PosA < PosC);
	asserta(PosC + CL <= L);

	uint Hi = PosC + CL - 1;
	uint PPL = Hi - PosA + 1;
	return PPL;
	}

void PDBChain::GetPalmCoreCoords(uint PosA, uint PosC, uint L,
  uint &Lo, uint &Hi)
	{
	asserta(PosA < PosC);
	asserta(PosC + CL <= L);

	Lo = (PosA <= PALMCOREM ? 0 : PosA - PALMCOREM);
	Hi = PosC + CL + PALMCOREM - 1;
	if (Hi >= L)
		Hi = L - 1;
	}

char PDBChain::GetMotifChar_Pos(uint Pos) const
	{
	if (SIZE(m_MotifPosVec) != 3)
		return '-';

	const uint PosA = m_MotifPosVec[A]; 
	const uint PosB = m_MotifPosVec[B]; 
	const uint PosC = m_MotifPosVec[C]; 

	if (Pos >= PosA && Pos < PosA + AL)
		return 'A';
	if (Pos >= PosB && Pos < PosB + BL)
		return 'B';
	if (Pos >= PosC && Pos < PosC + CL)
		return 'C';
	return '.';
	}

bool PDBChain::CheckMotifCoords(bool FailOnError) const
	{
	uint n = SIZE(m_MotifPosVec);
	if (n != 3)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), m_MotifPosVec.size()=%d", n);
		return false;
		}

	uint QL = SIZE(m_Seq);
	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	if (PosA + AL >= PosB)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), PosA+AL>PosB");
		return false;
		}

	if (PosB + BL >= PosC)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), PosB+BL>PosC");
		return false;
		}

	if (PosC + CL > QL)
		{
		if (FailOnError)
			Die("CheckMotifCoords(), PosC+CL > QC");
		return false;
		} 

	return true;
	}

bool PDBChain::CheckPPCMotifCoords(bool FailOnError) const
	{
	uint n = SIZE(m_MotifPosVec);
	if (n == 0)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), m_MotifPosVec empty");
		return false;
		}
	asserta(n == 3);

	uint QL = SIZE(m_Seq);
	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	if (PosA != 0)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), PosA=%u", PosA);
		return false;
		}

	if (PosA + AL >= PosB)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), PosA+AL>PosB");
		return false;
		}

	if (PosB + BL >= PosC)
		{
		if (FailOnError)
			Die("CheckPPCMotifCoords(), PosB+BL>PosC");
		return false;
		}

	if (PosC + CL != QL)
		{
		if (FailOnError)
			{
			LogMe();
			Die("CheckPPCMotifCoords(), PosC=%u +CL != QL=%u", PosC, QL);
			}
		return false;
		}

	return true;
	}

void PDBChain::GetPPC(PDBChain &PPC) const
	{
	PPC.Clear();

	asserta(SIZE(m_MotifPosVec) == 3);

	uint PosA = m_MotifPosVec[A];
	uint PosB = m_MotifPosVec[B];
	uint PosC = m_MotifPosVec[C];

	uint L = GetSeqLength();
	uint PPL = GetPalmPrintLength(PosA, PosC, L);

	asserta(PosB > PosA + AL);
	asserta(PosC > PosB + AL);
	asserta(PosC + CL <= L);

	uint PPC_PosA = 0;
	uint PPC_PosB = PosB - PosA;
	uint PPC_PosC = PosC - PosA;

	PPC.m_MotifPosVec.clear();
	PPC.m_MotifPosVec.push_back(PPC_PosA);
	PPC.m_MotifPosVec.push_back(PPC_PosB);
	PPC.m_MotifPosVec.push_back(PPC_PosC);

	vector<vector<double> > MotifCoords(3);
	GetPt(PosA, MotifCoords[A]);
	GetPt(PosB, MotifCoords[B]);
	GetPt(PosC, MotifCoords[C]);

	vector<double> t;
	vector<vector<double> > R;
	GetTriFormChain_tR(t, R, PPC);

	string MotifA, MotifB, MotifC;
	GetSubSeq(PosA, AL, MotifA);
	GetSubSeq(PosB, BL, MotifB);
	GetSubSeq(PosC, CL, MotifC);

	string MotifCoordStr;
	Ps(MotifCoordStr, "A:%u:%s B:%u:%s C:%u:%s",
	  PPC_PosA + 1, MotifA.c_str(),
	  PPC_PosB + 1, MotifB.c_str(),
	  PPC_PosC + 1, MotifC.c_str());

	PPC.m_Label = m_Label + " " + MotifCoordStr;
	PPC.CheckPPCMotifCoords(true);
	}

bool PDBChain::IsATOMLine(const string &Line)
	{
	if (Line[26] != ' ') // insertion code
		return false;
	if (strncmp(Line.c_str(), "ATOM  ", 6) == 0)
		return true;
	if (strncmp(Line.c_str(), "HETATM", 6) == 0)
		return true;
	return false;
	}

char PDBChain::GetChainCharFromATOMLine(const string &Line)
	{
	if (!IsATOMLine(Line))
		return 0;
	if (Line.size() < 22)
		return 0;
	return Line[21];
	}

double PDBChain::GetSmoothedCoord(uint Axis, uint i, uint N, uint w) const
	{
	double SumCoord = 0;
	double SumWeight = 0;
	double dN = (double) N;
	double dw = (double) w;
	int iQL = (int) GetSeqLength();
	double dQL = (double) iQL;
	for (int k = -1; k <= int(w); ++k)
		{
		double dPos = (int(i) + k)*dQL/dN;
		int iPos = int(dPos + 0.5);
		if (iPos < 0 || iPos >= iQL)
			continue;
		double Weight = dw - fabs(dPos - iPos);
		double Coord = GetCoord(Axis, (uint) iPos);
		SumCoord += Weight*Coord;
		SumWeight += Weight;
		}
	asserta(SumWeight > 0);
	double SmoothedCoord = SumCoord/SumWeight;
	return SmoothedCoord;
	}

char PDBChain::GetMotifA_D() const
	{
	asserta(SIZE(m_MotifPosVec) == 3);
	char D1 = m_Seq[m_MotifPosVec[A] + 3];
	if (D1 == 'D')
		return D1;
	char c1 = m_Seq[m_MotifPosVec[A] + 2];
	if (c1 == 'D')
		return 'D';
	char c2 = m_Seq[m_MotifPosVec[A] + 4];
	if (c2 == 'D')
		return 'D';
	return D1;
	}

char PDBChain::GetMotifB_G() const
	{
	asserta(SIZE(m_MotifPosVec) == 3);
	char G = m_Seq[m_MotifPosVec[B] + 0];
	if (G == 'G')
		return G;
	char c1 = m_Seq[m_MotifPosVec[B] + 1];
	if (c1 == 'G')
		return 'G';
	char c2 = m_Seq[m_MotifPosVec[B] + 2];
	if (c2 == 'G')
		return 'G';
	return G;
	}

char PDBChain::GetMotifC_D() const
	{
	asserta(SIZE(m_MotifPosVec) == 3);
	char D = m_Seq[m_MotifPosVec[C] + 2];
	if (D == 'D')
		return D;
	char c1 = m_Seq[m_MotifPosVec[C] + 3];
	if (c1 == 'D')
		return 'D';
	char c2 = m_Seq[m_MotifPosVec[C] + 4];
	if (c2 == 'D')
		return 'D';
	return D;
	}

void PDBChain::GetSuperMotif(string &s) const
	{
	s.clear();
	char D1 = GetMotifA_D();
	char G = GetMotifB_G();
	char D2 = GetMotifC_D();

	s += D1;
	s += G;
	s += D2;
	}

static uint GetDiffs3(const char *s, const char *t)
	{
	uint n = 0;
	for (uint i = 0; i < 3; ++i)
		if (s[i] != t[i])
			++n;
	return n;
	}

void PDBChain::GetBestMatchGDD(string &BestMatch) const
	{
	static const char *Motifs[] =
		{
		"GDD",
		"SDD",
		"ADD",
		"IDD",
		"GDN",
		};
	static const uint N = sizeof(Motifs)/sizeof(Motifs[0]);

	string CSeq;
	GetMotifSeq(C, CSeq);
	asserta(SIZE(CSeq) == 8);
	BestMatch = "XXX";
	uint BestDiffs = 3;
	for (uint i = 0; i < 5; ++i)
		{
		for (int j = 0; j < N; ++j)
			{
			uint Diffs = GetDiffs3(CSeq.c_str() + i, Motifs[j]);
			if (Diffs < BestDiffs)
				{
				BestMatch = CSeq.substr(i, 3);
				BestDiffs = Diffs;
				}
			}
		}
	}

const char *PDBChain::GetAcc(string &Acc) const
	{
	size_t n = m_Label.find(' ');
	if (n > 0)
		Acc = m_Label.substr(0, n);
	else
		Acc = m_Label;
	return Acc.c_str();
	}

void PDBChain::TruncateChain(uint Lo, uint Hi, PDBChain &Chain) const
	{
	Chain.Clear();
	asserta(Lo <= Hi);
	uint L = GetSeqLength();
	asserta(Hi < L);

	Chain.m_Label = m_Label;

	for (uint Pos = Lo; Pos <= Hi; ++Pos)
		{
#define c(x)	Chain.m_##x.push_back(m_##x[Pos])
		c(Seq);
		if (!m_ATOMs.empty())
			c(ATOMs);
		c(Xs);
		c(Ys);
		c(Zs);
#undef c
		}

	if (!m_MotifPosVec.empty())
		{
		asserta(SIZE(m_MotifPosVec) == 3);
		for (uint i = 0; i < 3; ++i)
			{
			uint Pos = m_MotifPosVec[i];
			asserta(Pos >= Lo);
			Chain.m_MotifPosVec.push_back(Pos - Lo);
			}
		asserta(SIZE(Chain.m_MotifPosVec) == 3);
		}
	}

void PDBChain::GetSphere(uint Pos, double Radius,
  uint MinPos, uint MaxPos, vector<uint> &PosVec) const
	{
	PosVec.clear();
	const uint L = GetSeqLength();
	asserta(MinPos <= MaxPos && MaxPos < L);
	for (uint Pos2 = MinPos; Pos2 <= MaxPos; ++Pos2)
		{
		if (Pos2 == Pos)
			continue;
		double d = GetDist(Pos, Pos2);
		if (d <= Radius)
			PosVec.push_back(Pos2);
		}
	}
