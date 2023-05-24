#ifndef pdbchain_h
#define pdbchain_h

#include "myutils.h"

class PDBChain
	{
public:
	string m_Label;
	string m_Seq;
	vector<double> m_Xs;
	vector<double> m_Ys;
	vector<double> m_Zs;
	vector<int> m_ResNrs;
	vector<vector<string> > m_ATOMs;
	vector<uint> m_MotifPosVec;
	string m_SS;

public:
	void Clear()
		{
		m_Label.clear();
		m_Seq.clear();
		m_Xs.clear();
		m_Ys.clear();
		m_Zs.clear();
		m_ResNrs.clear();
		m_MotifPosVec.clear();
		m_ATOMs.clear();
		m_SS.clear();
		}

	void SetMotifPosVec(uint PosA, uint PosB, uint PosC);
	uint GetSeqLength() const;
	void PrintSeqCoords(FILE *f) const;
	char FromPDBLines(const string &Label,
	  const vector<string> &Lines);
	void FromCal(const string &FileName);
	void FromCalLines(const vector<string> &Lines);
	void ParseCalLabelLine(const string &Line);
	void RenumberResidues(uint Start);
	void ToCal(FILE *f) const;
	void ToCalSeg(FILE *f, uint Pos, uint n) const;
	void ToCal(const string &FileName) const;
	void ToPDB(const string &FileName) const;
	void GetTriFormChain_tR(
	  const vector<double> &t,
	  const vector<vector<double> > &R,
	  PDBChain &XChain) const;
	void GetTriFormChain_MotifCoords(PDBChain &XChain) const;
	void GetTriFormChain_DGD(PDBChain &XChain) const;
	void LogMe(bool WithCoords = false) const;
	void GetMotifSeq(uint MotifIndex, string &s) const;
	void GetSubSeq(uint Pos, uint n, string &s) const;
	void GetXYZ(uint Pos, double &x, double &y, double &z) const;
	bool Get_CB_XYZ(uint Pos, double &x, double &y, double &z) const;
	bool Get_CB_Pt(uint Pos, vector<double> &Pt) const;
	double GetCoord(uint Axis, uint Pos) const;
	void GetPt(uint Pos, vector<double> &Pt) const;
	void SetPt(uint Pos, const vector<double> &Pt);
	double GetDist(uint Pos1, uint Pos2) const;
	double GetDist2(uint Pos1, uint Pos2) const;
	void GetMotifCoords(vector<vector<double> > &MotifCoords) const;
	void GetDGDCoords(vector<vector<double> > &DGDCoords) const;
	void GetMotifDists(double &AB, double &BC, double &AC) const;
	void GetMotifDists2(double &AB, double &BC, double &AC) const;
	uint GetMotifPos(uint MotifIndex) const;
	void GetSubSeq(uint StartPos, uint n,
	  bool FailOnOverflow, string &MotifSeq) const;
	void GetSS(string &SS) const;
	void SetSS()
		{
		GetSS(m_SS);
		}
	void GetPPC(PDBChain &PPC) const;
	void GetPC(PDBChain &PC) const;
	char GetMotifChar_Pos(uint Pos) const;
	bool CheckMotifCoords(bool FailOnError = true) const;
	bool CheckPPCMotifCoords(bool FailOnError = true) const;
	double GetSmoothedCoord(uint Axis, uint i, uint N, uint w) const;
	void GetSuperMotif(string &s) const;
	void GetBestMatchGDD(string &s) const;
	char GetMotifA_D() const;
	char GetMotifB_G() const;
	char GetMotifC_D() const;
	const char *GetAcc(string &Acc) const;
	void GetDistMx(uint Pos, uint L, vector<vector<double> > &Mx) const;
	void TruncateChain(uint Lo, uint Hi, PDBChain &Chain) const;
	void ToPML(FILE *f, const string &PDBFileName) const;
	void GetCAAtomLine(uint Pos, string &Line) const;
	int GetResidueNr(uint Pos) const;
	void GetResidueRange(uint PosLo, uint ResidueCount, int &ResLo,
	  int &ResHi) const;
	void GetSphere(uint Pos, double Radius,
	  uint MinPos, uint MaxPos,
	  vector<uint> &PosVec) const;

public:
	static uint GetMotifLength(uint MotifIndex);
	static uint GetPalmPrintLength(uint PosA, uint PosC, uint L);
	static void GetPalmCoreCoords(uint PosA, uint PosC, uint L,
	  uint &Lo, uint &Hi);
	static void ChainsFromLines(const string &Label,
	  const vector<string> &Lines, vector<PDBChain *> &Chains);
	static void ReadChainsFromFile(const string &FileName,
	  vector<PDBChain *> &Chains);
	static void AppendChainToLabel(string &Label, char Chain);
	static char GetChainCharFromATOMLine(const string &Line);
	static bool IsATOMLine(const string &Line);
	static int GetResidueNrFromATOMLine(const string &Line);
	static void HackHETAMLine(string &Line);
	static void GetAtomNameFromATOMLine(const string &Line,
	  string &AtomName);
	static void SetResidueNrInATOMLine(const string &InputLine,
	  uint ResidueNr, string &OutputLine);
	static void GetXYZFromATOMLine(const string &InputLine,
	  double &x, double &y, double &z);
	static void SetXYZInATOMLine(const string &InputLine,
	  double x, double y, double z, string &OutputLine);
	static bool GetFieldsFromResidueATOMLines(const vector<string> &Lines,
	  double &X, double &Y, double &Z, char &aa, int &ResNr);
	};

void ReadLinesFromFile(const string &FileName, vector<string> &Lines);
void ReadChains(const string &FileName,
  vector<PDBChain *> &Structures);
void ReadChains(const vector<string> &FileNames,
  vector<PDBChain *> &Structures);
void GetLabelFromFileName(const string &FileName, string &Label);

#endif // pdbchain_h
