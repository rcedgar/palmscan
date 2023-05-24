#pragma once

#include "trisearcher.h"
#include "pdbchain.h"

class TSHit
	{
public:
	PDBChain *m_Query = 0;
	const PDBChain *m_Ref = 0;
	uint m_QPosA = UINT_MAX;
	uint m_QPosB = UINT_MAX;
	uint m_QPosC = UINT_MAX;

	double m_TriRMSD2 = DBL_MAX;
	double m_MotifRMSD2 = DBL_MAX;

	double m_SketchPct = 0;
	string m_QSketch;
	string m_RSketch;
	string m_AnnotSketch;
	double m_Score = 0;

public:
	void Clear()
		{
		m_Query = 0;
		m_Ref = 0;
		m_QPosA = UINT_MAX;
		m_QPosB = UINT_MAX;
		m_QPosC = UINT_MAX;

		m_TriRMSD2 = DBL_MAX;
		m_MotifRMSD2 = DBL_MAX;

		m_SketchPct = 0;
		m_QSketch.clear();
		m_RSketch.clear();
		m_AnnotSketch.clear();

		m_Score = 0;
		}

	void WriteAln(FILE *f) const;
	void WriteTsv(FILE *f) const;
	void WriteSketch(FILE* f) const;
	void WriteReport(FILE *f) const;
	void WritePalmprintFasta(FILE *f) const;
	void SetSketch();
	};
