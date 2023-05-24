#pragma once

#include "pdbchain.h"

class GSProf
	{
public:
	string m_Label;
	string m_Seq;
	string m_Annot;
	vector<string> m_Names;
	vector<vector<double> > m_FeatureVec;

public:
	void Clear()
		{
		m_Label.clear();
		m_Seq.clear();
		m_Annot.clear();
		m_Names.clear();
		m_FeatureVec.clear();
		}

	void Init(const string &Label, const string &Seq,
	  const string &Annot);
	void AppendFeature(const string &Name,
	  const vector<double> &v);
	void ToTsv(FILE *f) const;
	bool FromTsv(FILE *f);
	bool FromTsv(const string &FileName);
	uint GetFeatureCount() const;
	uint GetSeqLength() const;
	};
