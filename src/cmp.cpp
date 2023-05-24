#include "myutils.h"
#include "cmp.h"
#include "abcxyz.h"

uint CMP::GetSeqPos(uint i, uint APos, uint BPos, uint CPos)
	{
	if (i >= CIX)
		return CPos + i - (CIX);
	if (i >= AL)
		return BPos + i - AL;
	return APos + i;
	}

char CMP::GetMotifChar(uint Ix)
	{
	if (Ix < AL)
		return 'A';
	if (Ix < AL + BL)
		return 'B';
	if (Ix < PPSPL)
		return 'C';
	asserta(false);
	return '!';
	}

void CMP::GetDistMx(const PDBChain &Q,
  uint APos, uint BPos, uint CPos,
  vector<vector<double> > &DistMx)
	{
	asserta(APos != UINT_MAX);
	asserta(BPos != UINT_MAX);
	asserta(CPos != UINT_MAX);

	DistMx.clear();
	const uint N = CIX + CL;
	DistMx.resize(N);
	for (uint i = 0; i < N; ++i)
		DistMx[i].resize(N, DBL_MAX);

	for (uint i = 0; i < N; ++i)
		{
		uint SeqPosi = GetSeqPos(i, APos, BPos, CPos);
		for (uint j = 0; j < N; ++j)
			{
			uint SeqPosj = GetSeqPos(j, APos, BPos, CPos);
			double d = Q.GetDist(SeqPosi, SeqPosj);
			DistMx[i][j] = d;
			}
		}
	}

void CMP::MxToFile(FILE *f, const string &Name,
  const vector<vector<double> > &Mx) const
	{
	if (f == 0)
		return;

	asserta(SIZE(Mx) == PPSPL);
	for (uint i = 1; i < PPSPL; ++i)
		{
		fprintf(f, "%s\t%u", Name.c_str(), i);
		for (uint j = 0; j <= i; ++j)
			{
			double StdDev = Mx[i][j];
			if (i == j)
				asserta(StdDev == 0);
			asserta(Mx[j][i] == Mx[i][j]);
			fprintf(f, "\t%.3g", StdDev);
			}
		fprintf(f, "\n");
		}
	}

void CMP::ToFile(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	asserta(f != 0);
	const uint N = SIZE(m_DistMxVec);
	asserta(SIZE(m_RefLabels) == N);
	fprintf(f, "CMP\t%u", N);
	if (m_Label == "")
		fprintf(f, "\t.\n");
	else
		fprintf(f, "\t%s\n", m_Label.c_str());
	MxToFile(f, "mean", m_RefMeans);
	MxToFile(f, "stddev", m_StdDevs);
	for (uint i = 0; i < N; ++i)
		{
		string Name = "ref." + m_RefLabels[i];
		MxToFile(f, Name, m_DistMxVec[i]);
		}
	CloseStdioFile(f);
	}

void CMP::MxFromFile(FILE *f, string &Name,
  vector<vector<double> > &Mx)
	{
	Mx.resize(PPSPL);
	for (uint i = 0; i < PPSPL; ++i)
		Mx[i].resize(PPSPL, DBL_MAX);

	string Line;
	vector<string> Fields;
	for (uint i = 1; i < PPSPL; ++i)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			Die("Premature EOF in CMP file");
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == i+3);
		if (i == 1)
			Name = Fields[0];
		else
			asserta(Fields[0] == Name);
		asserta(StrToUint(Fields[1]) == i);

		for (uint j = 0; j <= i; ++j)
			{
			double Value = StrToFloat(Fields[j+2]);
			if (i == j)
				asserta(Value == 0);
			Mx[i][j] = Value;
			Mx[j][i] = Value;
			}
		}
	}

void CMP::FromFile(FILE *f)
	{
	Clear();
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	if (!Ok)
		Die("Premature EOF in CMP file");
	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 3 || Fields[0] != "CMP")
		Die("Invalid CMP file (hdr)");
	const uint N = StrToUint(Fields[1]);
	m_Label = Fields[2];

	string Name;
	MxFromFile(f, Name, m_RefMeans);
	asserta(Name == "mean");

	MxFromFile(f, Name, m_StdDevs);
	asserta(Name == "stddev");

	m_DistMxVec.resize(N);
	for (uint i = 0; i < N; ++i)
		{
		MxFromFile(f, Name, m_DistMxVec[i]);
		asserta(Name.substr(0, 4) == "ref.");
		m_RefLabels.push_back(Name.substr(4, string::npos));
		}
	}

void CMP::FromFile(const string &FileName)
	{
	asserta(FileName != "");
	FILE *f = OpenStdioFile(FileName);
	FromFile(f);
	CloseStdioFile(f);
	}
