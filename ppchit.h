#pragma once

#include "trisearcher.h"
#include "pdbchain.h"

class PpcHit
	{
public:
	const PDBChain *m_Query = 0;
	const PDBChain *m_Ref = 0;
	double m_MotifRMSD2 = DBL_MAX;
	double m_SketchPct = 0;
	string m_QSketch;
	string m_RSketch;
	string m_AnnotSketch;

public:
	void Clear()
		{
		m_Query = 0;
		m_Ref = 0;
		m_MotifRMSD2 = DBL_MAX;
		m_SketchPct = 0;
		m_QSketch.clear();
		m_RSketch.clear();
		m_AnnotSketch.clear();
		}

	void WriteAln(FILE *f) const;
	void WriteTsv(FILE *f) const;
	void WritePalmprintFasta(FILE *f) const;
	void WritePalmprintPDB(const string &FileNamePrefix) const;
	void SetSketch();
	void WriteSketch(FILE* f) const;
	double GetScore() const;
	};
