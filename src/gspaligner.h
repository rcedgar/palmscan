#pragma once

#include "gsprof.h"

class GSPAligner
	{
public:
	const GSProf *m_Q = 0;
	const GSProf *m_T = 0;
	uint m_FeatureCount = 0;
	double m_Bias = 0.2;
	uint m_MxSize = 0;
	float **m_Mx = 0;

public:
	void Clear()
		{
		m_Q = 0;
		m_T = 0;
		}

	void Align(const GSProf &Q, const GSProf &T);
	double GetScore_Pos(uint QPos, uint TPos) const;
	double GetScore_Vecs(const vector<double>  &QVec,
	  const vector<double>  &TVec) const;
	void Alloc(uint QL, uint TL);
	void SetMx();
	void MxToTsv(FILE *f) const;
	};
