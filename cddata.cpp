#include "myutils.h"
#include "cddata.h"

void CDData::ToTsv(const string &Name, FILE *f) const
	{
	if (f == 0)
		return;

	asserta(m_Info != 0);
	const uint Size = m_Info->GetSize();
	fprintf(f, "cddata\t%s\t%u\n", Name.c_str(), Size);
	for (uint Ix1 = 0; Ix1 < Size; ++Ix1)
		{
		uint MotifIndex1, i1;
		m_Info->GetCoord(Ix1, MotifIndex1, i1);
		const string &MotifName1 =
		  m_Info->m_MotifNames[MotifIndex1];
		for (uint Ix2 = 0; Ix2 < Ix1; ++Ix2)
			{
			uint MotifIndex2, i2;
			m_Info->GetCoord(Ix2, MotifIndex2, i2);
			double Value = GetByIx(Ix1, Ix2);
			const string &MotifName2 =
			  m_Info->m_MotifNames[MotifIndex2];
			fprintf(f, "%s\t%u\t%s\t%u\t%.4g\n",
			  MotifName1.c_str(), i1,
			  MotifName2.c_str(), i2,
			  Value);
			}
		}
	}

void CDData::FromTsv(const CDInfo &Info, FILE *f)
	{
	Init(Info);

	string Line;
	vector<string> Fields;

	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 3);
	asserta(Fields[0] == "cddata");

	const uint Size = m_Info->GetSize();
	for (uint Ix1 = 0; Ix1 < Size; ++Ix1)
		{
		uint MotifIndex1, i1;
		m_Info->GetCoord(Ix1, MotifIndex1, i1);
		const string &MotifName1 =
		  m_Info->m_MotifNames[MotifIndex1];
		for (uint Ix2 = 0; Ix2 < Ix1; ++Ix2)
			{
			uint MotifIndex2, i2;
			m_Info->GetCoord(Ix2, MotifIndex2, i2);
			const string &MotifName2 =
			  m_Info->m_MotifNames[MotifIndex2];
			bool Ok = ReadLineStdioFile(f, Line);
			asserta(Ok);
			Split(Line, Fields, '\t');
			asserta(SIZE(Fields) == 5);

			asserta(Fields[0] == MotifName1);
			asserta(StrToUint(Fields[1]) == i1);
			asserta(Fields[2] == MotifName2);
			asserta(StrToUint(Fields[3]) == i2);

			double Value = StrToFloat(Fields[4]);
			SetByIx(Ix1, Ix2, Value);
			}
		}
	}

void CDData::Init(const CDInfo &Info)
	{
	m_Info = &Info;
	uint Size = m_Info->GetSize();
	m_Data.clear();
	m_Data.resize(Size);
	for (uint i = 0; i < Size; ++i)
		m_Data[i].resize(Size, DBL_MAX);
	}

double CDData::Get(uint MotifIndex1, uint i1,
  uint MotifIndex2, uint i2) const
	{
	uint Ix1 = m_Info->GetIx(MotifIndex1, i1);
	uint Ix2 = m_Info->GetIx(MotifIndex2, i2);
	double Value = GetByIx(Ix1, Ix2);
	return Value;
	}

void CDData::Set(uint MotifIndex1, uint i1,
  uint MotifIndex2, uint i2, double Value)
	{
	uint Ix1 = m_Info->GetIx(MotifIndex1, i1);
	uint Ix2 = m_Info->GetIx(MotifIndex2, i2);
	SetByIx(Ix1, Ix2, Value);
	}

double CDData::GetByIx(uint Ix1, uint Ix2) const
	{
	assert(Ix1 < SIZE(m_Data));
	assert(Ix2 < SIZE(m_Data[Ix1]));
	double Value = m_Data[Ix1][Ix2];
	return Value;
	}

void CDData::SetByIx(uint Ix1, uint Ix2, double Value)
	{
	assert(Ix1 < SIZE(m_Data));
	assert(Ix2 < SIZE(m_Data[Ix1]));
	m_Data[Ix1][Ix2] = Value;
	m_Data[Ix2][Ix1] = Value;
	}
