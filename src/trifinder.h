#pragma once

#include "pdbchain.h"

class TriFinder
	{
public:

	const PDBChain *m_Query = 0;
	uint m_PosA = UINT_MAX;
	uint m_PosB = UINT_MAX;
	uint m_PosC = UINT_MAX;
	double m_Score = 0;
	uint m_TriCount = 0;

public:
	static uint L_ABmin;
	static uint L_ACmin;
	static uint L_BCmin;

	static uint L_ABmax;
	static uint L_ACmax;
	static uint L_BCmax;

	static double d_ABmin;
	static double d_ACmin;
	static double d_BCmin;

	static double d_ABmax;
	static double d_ACmax;
	static double d_BCmax;

public:
	void Clear()
		{
		m_Query = 0;
		m_PosA = UINT_MAX;
		m_PosB = UINT_MAX;
		m_PosC = UINT_MAX;
		m_Score = 0;
		m_TriCount = 0;
		}

	void LogMe(FILE *f = g_fLog) const;
	void Find(const PDBChain &Query, bool DGD);
	void LogTri(uint PosA, uint PosB, uint PosC) const;

public:
	static void SetSigmas(double Sigmas);
	};
