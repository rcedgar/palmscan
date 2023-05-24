#pragma once

#include "cmpsearcher.h"
#include "pdbchain.h"
#include "cmp.h"

class DSHit;

class CMSearcher
	{
public:
	CMP m_ProfDB;

public:
	CMPSearcher m_PS;

	PDBChain *m_Query = 0;

	uint m_RefIndex = UINT_MAX;
	uint m_PosA = UINT_MAX;
	uint m_PosB = UINT_MAX;
	uint m_PosC = UINT_MAX;
	double m_Score = 0;

public:
	void ClearSearch()
		{
		m_Query = 0;
		m_RefIndex = UINT_MAX;
		m_PosA = UINT_MAX;
		m_PosB = UINT_MAX;
		m_PosC = UINT_MAX;
		m_Score = DBL_MAX;
		}

	void Search(PDBChain &Query);
	uint GetPSSMStarts(uint &PosA, uint &PosB, uint &PosC,
	  double &Score) const;
	};
