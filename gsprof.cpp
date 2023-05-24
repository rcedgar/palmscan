#include "myutils.h"
#include "gsprof.h"

void GSProf::Init(const string &Label, const string &Seq,
  const string &Annot)
	{
	Clear();
	m_Label = Label;
	m_Seq = Seq;
	m_Annot = Annot;
	}

void GSProf::AppendFeature(const string &Name,
  const vector<double> &v)
	{
	m_Names.push_back(Name);
	const uint L = GetSeqLength();
	asserta(SIZE(v) == L);
	if (m_FeatureVec.empty())
		m_FeatureVec.resize(L);
	for (uint i = 0; i < L; ++i)
		m_FeatureVec[i].push_back(v[i]);
	}

uint GSProf::GetFeatureCount() const
	{
	const uint N = SIZE(m_Names);
	return N;
	}

uint GSProf::GetSeqLength() const
	{
	uint L = SIZE(m_Seq);
	asserta(SIZE(m_Annot) == L);
	return L;
	}

void GSProf::ToTsv(FILE *f) const
	{
	const uint NF = GetFeatureCount();
	const uint L = GetSeqLength();
	const char *Lab = m_Label.c_str();

	fprintf(f, "%s", Lab);
	fprintf(f, "\t%u", L);
	fprintf(f, "\tSeq");
	fprintf(f, "\tAnnot");
	for (uint i = 0; i < NF; ++i)
		fprintf(f, "\t%s", m_Names[i].c_str());
	fprintf(f, "\n");

	for (uint Pos = 0; Pos < L; ++Pos)
		{
		fprintf(f, "%s", Lab);
		fprintf(f, "\t%u", Pos);
		fprintf(f, "\t%c", m_Seq[Pos]);
		fprintf(f, "\t%c", m_Annot[Pos]);
		for (uint i = 0; i < NF; ++i)
			fprintf(f, "\t%.3g", m_FeatureVec[Pos][i]);
		fprintf(f, "\n");
		}
	}

bool GSProf::FromTsv(const string &FileName)
	{
	if (FileName == "")
		return false;
	FILE *f = OpenStdioFile(FileName);
	bool Ok = FromTsv(f);
	CloseStdioFile(f);
	return Ok;
	}

bool GSProf::FromTsv(FILE *f)
	{
	Clear();
	string Line;
	bool Ok = ReadLineStdioFile(f, Line);
	if (!Ok)
		return false;

	vector<string> Fields;
	Split(Line, Fields, '\t');
	const uint N = SIZE(Fields);
	asserta(N > 2);
	m_Label = Fields[0];
	const uint L = StrToUint(Fields[1]);
	asserta(Fields[2] == "Seq");
	asserta(Fields[3] == "Annot");
	for (uint i = 4; i < N; ++i)
		m_Names.push_back(Fields[i]);

	uint NF = SIZE(m_Names);
	m_FeatureVec.resize(L);
	for (uint Pos = 0; Pos < L; ++Pos)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == N);
		asserta(Fields[0] == m_Label);
		asserta(StrToUint(Fields[1]) == Pos);
		asserta(SIZE(Fields[2]) == 1);
		asserta(SIZE(Fields[3]) == 1);
		char c = Fields[2][0];
		char a = Fields[3][0];
		m_Seq += c;
		m_Annot += a;

		for (uint i = 0; i < NF; ++i)
			{
			double x = StrToFloat(Fields[4+i]);
			m_FeatureVec[Pos].push_back(x);
			}
		}

	return true;
	}
