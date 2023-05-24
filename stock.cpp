#include "myutils.h"
#include "stock.h"

void Stock::FromFile(const string &FileName)
	{
	FILE *f = OpenStdioFile(FileName);

	uint BlockIndex = 0;
	uint SeqIndex = 0;
	string Line;
	vector<string> Fields;
	for (;;)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		asserta(Ok);
		if (Line == "//")
			break;
		if (Line == "")
			{
			++BlockIndex;
			SeqIndex = 0;
			continue;
			}

		SplitWhite(Line, Fields);
		uint FieldCount = SIZE(Fields);
		asserta(FieldCount > 1);
		if (Line[0] == '#')
			{
		// Annotation line
			if (BlockIndex == 0)
				m_Header.push_back(Line);
			else
				{
				asserta(FieldCount == 3);
				const string &Feature = Fields[1];
				const string &Row = Fields[2];
				if (BlockIndex == 1)
					{
					map<string, string>::iterator p =
					  m_FeatureToRow.find(Feature);
					asserta(p == m_FeatureToRow.end());
					m_FeatureToRow[Feature] = Row;
					}
				else
					{
					map<string, string>::iterator p =
					  m_FeatureToRow.find(Feature);
					asserta(p != m_FeatureToRow.end());
					p->second += Row;
					}
				}
			}
		else
			{
		// Sequence line
			asserta(FieldCount == 2);
			const string &Label = Fields[0];
			const string &Row = Fields[1];
			asserta(BlockIndex != 0);
			if (BlockIndex == 1)
				{
				m_Labels.push_back(Label);
				m_Rows.push_back(Row);
				}
			else
				{
				asserta(SeqIndex < SIZE(m_Labels));
				asserta(m_Labels[SeqIndex] == Label);
				m_Rows[SeqIndex] += Row;
				}
			++SeqIndex;
			}
		}
	CloseStdioFile(f);

	const uint SeqCount = SIZE(m_Labels);
	asserta(SIZE(m_Rows) == SeqCount);
	asserta(SeqCount > 0);
	const uint ColCount = SIZE(m_Rows[0]);
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		asserta(SIZE(m_Rows[SeqIndex]) == ColCount);

	for (map<string, string>::const_iterator p = m_FeatureToRow.begin();
	  p != m_FeatureToRow.end(); ++p)
		{
		const string &Row = p->second;
		asserta(SIZE(Row) == ColCount);
		}
	}

void Stock::ToFasta(FILE *f) const
	{
	if (f == 0)
		return;

	const uint SeqCount = GetSeqCount();
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Label = m_Labels[SeqIndex];
		const string &Seq = m_Rows[SeqIndex];
		SeqToFasta(f, Label, Seq);
		}
	}

void Stock::GetMotifCols(uint &ColA, uint &ColB, uint &ColC) const
	{
	ColA = UINT_MAX;
	ColB = UINT_MAX;
	ColC = UINT_MAX;

	uint NA = 0;
	uint NB = 0;
	uint NC = 0;
	map<string, string>::const_iterator p = m_FeatureToRow.find("ABC");
	asserta(p != m_FeatureToRow.end());
	const string &Row = p->second;
	
	const uint ColCount = SIZE(Row);
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		switch (Row[Col])
			{
		case 'A':
			if (NA == 0)
				ColA = Col;
			asserta(Col == ColA + NA);
			++NA;
			break;
		case 'B':
			if (NB == 0)
				ColB = Col;
			asserta(Col == ColB + NB);
			++NB;
			break;
		case 'C':
			if (NC == 0)
				ColC = Col;
			asserta(Col == ColC + NC);
			++NC;
			break;
			}
		}
	asserta(NA == 12);
	asserta(NB == 14);
	asserta(NC == 8);
	}

void Stock::GetMotifsStrings(vector<string> &As, vector<string> &Bs,
  vector<string> &Cs) const
	{
	As.clear();
	Bs.clear();
	Cs.clear();

	uint ColA, ColB, ColC;
	GetMotifCols(ColA, ColB, ColC);

	const uint SeqCount = SIZE(m_Rows);
	uint ColCount = UINT_MAX;
	for (uint SeqIndex = 0; SeqIndex < SeqCount; ++SeqIndex)
		{
		const string &Row = m_Rows[SeqIndex];
		if (SeqIndex == 0)
			ColCount = SIZE(Row);
		else
			asserta(SIZE(Row) == ColCount);

		asserta(ColA + 12 <= ColCount);
		asserta(ColB + 14 <= ColCount);
		asserta(ColC + 8 <= ColCount);

		string A, B, C;
		for (uint i = 0; i < 12; ++i)
			A += Row[ColA + i];
		for (uint i = 0; i < 14; ++i)
			B += Row[ColB + i];
		for (uint i = 0; i < 8; ++i)
			C += Row[ColC + i];

		As.push_back(A);
		Bs.push_back(B);
		Cs.push_back(C);
		}
	}
