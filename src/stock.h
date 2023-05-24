#pragma once

#include <map>

class Stock
	{
public:
	vector<string> m_Header;
	vector<string> m_Labels;
	vector<string> m_Rows;
	map<string, string> m_FeatureToRow;

public:
	void Clear()
		{
		m_Header.clear();
		m_Labels.clear();
		m_Rows.clear();
		m_FeatureToRow.clear();
		}

	uint GetSeqCount() const { return SIZE(m_Labels); }
	void FromFile(const string &FileName);
	void ToFasta(FILE *f) const;
	void GetMotifCols(uint &ColA, uint &ColB, uint &ColC) const;
	void GetMotifsStrings(vector<string> &As, vector<string> &Bs,
	  vector<string> &Cs) const;
	};
