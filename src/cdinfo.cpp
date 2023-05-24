#include "myutils.h"
#include "cdinfo.h"

void CDInfo::LogMe() const
	{
	Log("%u motifs: \n", m_MotifCount);
	for (uint i = 0; i < m_MotifCount; ++i)
		Log(" %s(%u)",
		  m_MotifNames[i].c_str(),
		  m_MotifLengths[i]);
	}

void CDInfo::Init(const vector<string> &MotifNames,
	vector<uint> &MotifLengths)
	{
	Clear();
	m_MotifCount = SIZE(MotifNames);
	asserta(SIZE(MotifLengths) == m_MotifCount);
	m_MotifNames = MotifNames;
	m_MotifLengths = MotifLengths;
	uint Offset = 0;
	for (uint i = 0; i < m_MotifCount; ++i)
		{
		uint ML = m_MotifLengths[i];
		m_Offsets.push_back(Offset);
		Offset += ML;
		}
	m_Offsets.push_back(Offset);

	uint Size = GetSize();
	for (uint Ix = 0; Ix < Size; ++Ix)
		{
		uint MotifIndex, i;
		GetCoord(Ix, MotifIndex, i);
		m_MotifIndexes.push_back(MotifIndex);
		m_MotifCoords.push_back(i);
		}
	}

void CDInfo::GetCoord(uint Ix, uint &MotifIndex, uint &i) const
	{
	for (MotifIndex = 0; ; ++MotifIndex)
		{
		if (Ix < m_Offsets[MotifIndex+1])
			{
			i = Ix - m_Offsets[MotifIndex];
			return;
			}
		}
	asserta(false);
	}

uint CDInfo::GetIx(uint MotifIndex, uint i) const
	{
	assert(MotifIndex < SIZE(m_Offsets));
	assert(i < m_MotifLengths[i]);
	uint Ix = m_Offsets[MotifIndex] + i;
	return Ix;
	}

void CDInfo::ToTsv(FILE *f) const
	{
	if (f == 0)
		return;
	fprintf(f, "cdinfo\t%u\n", m_MotifCount);
	for (uint i = 0; i < m_MotifCount; ++i)
		fprintf(f, "%u\t%u\t%s\n",
		  i,
		  m_MotifLengths[i],
		  m_MotifNames[i].c_str());
	}

void CDInfo::FromTsv(FILE *f)
	{
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	asserta(Ok);
	Split(Line, Fields, '\t');
	asserta(SIZE(Fields) == 2);
	asserta(Fields[0] == "cdinfo");
	uint MotifCount = StrToUint(Fields[1]);
	asserta(MotifCount > 0 && MotifCount < 99);
	vector<string> Names;
	vector<uint> Lengths;
	for (uint i = 0; i < MotifCount; ++i)
		{
		Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == 3);
		asserta(StrToUint(Fields[0]) == i);
		uint ML = StrToUint(Fields[1]);
		asserta(ML > 0 && ML < 99);
		Lengths.push_back(ML);
		Names.push_back(Fields[2]);
		}
	Init(Names, Lengths);
	}

const char *CDInfo::GetMotifName(uint MotifIndex) const
	{
	asserta(MotifIndex < SIZE(m_MotifNames));
	return m_MotifNames[MotifIndex].c_str();
	}

uint CDInfo::GetMotifLength(uint MotifIndex) const
	{
	asserta(MotifIndex < SIZE(m_MotifLengths));
	return m_MotifLengths[MotifIndex];
	}

uint CDInfo::GetSize() const
	{
	assert(m_MotifCount > 0 && SIZE(m_Offsets) == m_MotifCount+1);
	uint ML = m_MotifLengths[m_MotifCount-1];
	uint Off = m_Offsets[m_MotifCount-1];
	uint Size = ML + Off;
	return Size;
	}
